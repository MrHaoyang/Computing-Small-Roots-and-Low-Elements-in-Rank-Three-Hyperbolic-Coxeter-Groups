"""Algebraic source of truth for rank-three hyperbolic Coxeter systems."""

from collections import deque
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import combinations, product
from math import cos, inf, isinf, pi, sqrt
import re

import numpy as np


GENERATOR_ORDER = ("1", "2", "3")
INFINITE_ORDER = inf
INFINITE_ORDER_TOKENS = {"inf", "infinity", "∞", "unbounded"}


def normalize_coxeter_order(value):
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in INFINITE_ORDER_TOKENS:
            return INFINITE_ORDER
        try:
            value = int(normalized)
        except ValueError as exc:
            raise ValueError(f"Invalid Coxeter order: {value}") from exc
    if value == INFINITE_ORDER:
        return INFINITE_ORDER
    if not isinstance(value, int) or isinstance(value, bool) or value < 2:
        raise ValueError("A Coxeter order must be an integer at least 2 or infinity.")
    return value


def coxeter_order_reciprocal(value):
    value = normalize_coxeter_order(value)
    return Fraction(0, 1) if isinf(value) else Fraction(1, value)


def coxeter_order_cosine(value):
    value = normalize_coxeter_order(value)
    return 1.0 if isinf(value) else cos(pi / value)


def coxeter_order_label(value):
    value = normalize_coxeter_order(value)
    return "∞" if isinf(value) else str(value)


def serialized_coxeter_order(value):
    value = normalize_coxeter_order(value)
    return "inf" if isinf(value) else value


def lorentz_normals_from_gram(gram, tolerance=1e-9):
    """Factor a signature-(2,1) Gram matrix into Lorentzian normals."""
    gram = np.asarray(gram, dtype=np.float64)
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    negative = [index for index, value in enumerate(eigenvalues) if value < -tolerance]
    positive = [index for index, value in enumerate(eigenvalues) if value > tolerance]
    if len(negative) != 1 or len(positive) != 2:
        raise ValueError("The Gram matrix does not have Lorentz signature (2,1).")
    order = negative + positive
    ordered_values = eigenvalues[order]
    ordered_vectors = eigenvectors[:, order]
    normal_matrix = np.diag(np.sqrt(np.abs(ordered_values))) @ ordered_vectors.T

    minkowski = np.diag((-1.0, 1.0, 1.0))
    equations = normal_matrix.T @ minkowski
    center = np.linalg.solve(equations, -np.ones(3))
    norm = float(center @ minkowski @ center)
    if norm >= -tolerance:
        raise ArithmeticError("The center of the fundamental chamber is not on the hyperboloid.")
    if center[0] < 0:
        normal_matrix[0, :] *= -1.0

    np.testing.assert_allclose(
        normal_matrix.T @ minkowski @ normal_matrix,
        gram,
        atol=max(tolerance * 10, 1e-10),
    )
    return tuple(normal_matrix[:, index].copy() for index in range(3))


@dataclass
class RootRecord:
    key: tuple[int, ...]
    coefficients: np.ndarray
    reflection: np.ndarray
    depth: int
    is_small: bool = False


@dataclass
class WallRecord:
    root_key: tuple[int, ...]
    reflection: np.ndarray
    wall_type: str
    chamber_key: tuple[int, ...]
    neighbor_key: tuple[int, ...]


@dataclass
class ChamberRecord:
    key: tuple[int, ...]
    word: str
    matrix: np.ndarray
    aliases: set[str] = field(default_factory=set)
    length: int = 0
    inversions: tuple[tuple[int, ...], ...] = ()
    walls: dict[str, WallRecord] = field(default_factory=dict)


class CoxeterSystem3:
    """The geometric representation of a rank-three triangle Coxeter group.

    Roots use column coordinates in the ordered simple-root basis
    ``(alpha_1, alpha_2, alpha_3)``.  A displayed word ``123`` has matrix
    ``R_1 @ R_2 @ R_3``.
    """

    def __init__(self, m12, m13, m23, tolerance=1e-9):
        parameters = tuple(
            normalize_coxeter_order(value)
            for value in (m12, m13, m23)
        )
        if sum(coxeter_order_reciprocal(value) for value in parameters) >= 1:
            raise ValueError("The parameters do not satisfy the hyperbolicity condition: 1/m12 + 1/m13 + 1/m23 < 1.")

        self.m12, self.m13, self.m23 = parameters
        self.tolerance = float(tolerance)
        self.max_element_residual = 0.0
        self.max_root_residual = 0.0
        self._element_buckets = {}
        self._element_entries = []
        self._root_key_buckets = {}
        self._root_key_entries = []
        self.gram = np.array(
            [
                [
                    1.0,
                    -coxeter_order_cosine(self.m12),
                    -coxeter_order_cosine(self.m13),
                ],
                [
                    -coxeter_order_cosine(self.m12),
                    1.0,
                    -coxeter_order_cosine(self.m23),
                ],
                [
                    -coxeter_order_cosine(self.m13),
                    -coxeter_order_cosine(self.m23),
                    1.0,
                ],
            ],
            dtype=np.float64,
        )

        eigenvalues = np.linalg.eigvalsh(self.gram)
        positive = int(np.count_nonzero(eigenvalues > self.tolerance))
        negative = int(np.count_nonzero(eigenvalues < -self.tolerance))
        zero = 3 - positive - negative
        self.signature = (positive, negative, zero)
        if self.signature != (2, 1, 0):
            raise ValueError(
                "The Gram matrix does not have rank-three hyperbolic signature (2,1): "
                f"received {self.signature}."
            )

        identity = np.eye(3, dtype=np.float64)
        self.simple_roots = tuple(identity[:, index].copy() for index in range(3))
        self.generator_matrices = {}
        for index, generator in enumerate(GENERATOR_ORDER):
            simple_root = self.simple_roots[index]
            self.generator_matrices[generator] = identity - 2.0 * np.outer(
                simple_root,
                self.gram[index, :],
            )
        self._root_registry = {}

    @staticmethod
    def _normalize_word(word):
        normalized = str(word).strip().lower()
        if normalized == "e":
            return ""
        if not re.fullmatch(r"[123]*", normalized):
            raise ValueError(f"Invalid Coxeter word: {word}")
        return normalized

    def word_matrix(self, word):
        matrix = np.eye(3, dtype=np.float64)
        for generator in self._normalize_word(word):
            matrix = matrix @ self.generator_matrices[generator]
        return matrix

    def element_key(self, word_or_matrix):
        matrix = (
            self.word_matrix(word_or_matrix)
            if isinstance(word_or_matrix, str)
            else np.asarray(word_or_matrix, dtype=np.float64)
        )
        if matrix.shape != (3, 3):
            raise ValueError("A group element matrix must be 3x3.")
        coordinates = self._key_coordinates(matrix)
        quantum = max(self.tolerance * 1000, 1e-6)
        bucket = tuple(
            int(value)
            for value in np.rint(coordinates / quantum).astype(np.int64)
        )
        candidates = self._bucket_candidates(
            coordinates,
            quantum,
            self._element_buckets,
        )
        match = self._find_residual_match(matrix, candidates, is_root=False)
        if match is not None:
            primary = self._element_buckets.setdefault(bucket, [])
            if not any(existing is match for existing in primary):
                primary.append(match)
            return match[1]

        key = bucket + (len(self._element_entries),)
        entry = (matrix.copy(), key)
        self._element_buckets.setdefault(bucket, []).append(entry)
        self._element_entries.append(entry)
        return key

    @staticmethod
    def _key_coordinates(value):
        flattened = np.asarray(value, dtype=np.float64).reshape(-1)
        scale = max(1.0, float(np.max(np.abs(flattened))))
        return np.concatenate((flattened / scale, (np.log(scale),)))

    def _bucket_candidates(self, value, quantum, buckets):
        options = []
        threshold = max(self.tolerance * 20, 1e-8)
        for coordinate in np.asarray(value).reshape(-1):
            scaled = float(coordinate) / quantum
            primary = int(np.rint(scaled))
            coordinate_options = [primary]
            if (scaled - (primary - 0.5)) * quantum <= threshold:
                coordinate_options.append(primary - 1)
            if ((primary + 0.5) - scaled) * quantum <= threshold:
                coordinate_options.append(primary + 1)
            options.append(tuple(dict.fromkeys(coordinate_options)))

        candidates = []
        seen = set()
        for bucket in product(*options):
            for entry in buckets.get(tuple(bucket), ()):
                identity = id(entry)
                if identity not in seen:
                    seen.add(identity)
                    candidates.append(entry)
        return candidates

    def _find_residual_match(self, value, entries, is_root):
        if not entries:
            return None
        values = np.stack([entry[0] for entry in entries])
        absolute_residuals = np.max(
            np.abs(values - value),
            axis=tuple(range(1, values.ndim)),
        )
        entry_scales = np.max(
            np.abs(values),
            axis=tuple(range(1, values.ndim)),
        )
        value_scale = float(np.max(np.abs(value)))
        normalized_residuals = absolute_residuals / (
            1.0 + np.maximum(entry_scales, value_scale)
        )
        close = normalized_residuals <= self.tolerance * 10
        indexes = np.flatnonzero(close)
        if not indexes.size:
            return None
        match_index = int(indexes[np.argmin(normalized_residuals[indexes])])
        entry = entries[match_index]
        residual = float(normalized_residuals[match_index])
        if is_root:
            self.max_root_residual = max(self.max_root_residual, residual)
        else:
            self.max_element_residual = max(self.max_element_residual, residual)
        return entry

    def same_element(self, left, right):
        left_matrix = self.word_matrix(left)
        right_matrix = self.word_matrix(right)
        residual = np.max(np.abs(left_matrix - right_matrix)) / (
            1.0
            + max(
                np.max(np.abs(left_matrix)),
                np.max(np.abs(right_matrix)),
            )
        )
        return bool(residual <= self.tolerance * 10)

    def _positive_root(self, root):
        coefficients = np.asarray(root, dtype=np.float64).reshape(3).copy()
        coefficients[np.abs(coefficients) <= self.tolerance * 10] = 0.0
        if np.all(coefficients >= 0.0):
            return coefficients
        if np.all(coefficients <= 0.0):
            return -coefficients
        raise ArithmeticError(
            "The computed real root is neither positive nor negative; the generator order or numerical error is inconsistent."
        )

    def root_key(self, root):
        positive = self._positive_root(root)
        coordinates = self._key_coordinates(positive)
        quantum = max(self.tolerance * 1000, 1e-6)
        bucket = tuple(
            int(value)
            for value in np.rint(coordinates / quantum).astype(np.int64)
        )
        candidates = self._bucket_candidates(
            coordinates,
            quantum,
            self._root_key_buckets,
        )
        match = self._find_residual_match(positive, candidates, is_root=True)
        if match is not None:
            primary = self._root_key_buckets.setdefault(bucket, [])
            if not any(existing is match for existing in primary):
                primary.append(match)
            return match[1]

        key = bucket + (len(self._root_key_entries),)
        entry = (positive.copy(), key)
        self._root_key_buckets.setdefault(bucket, []).append(entry)
        self._root_key_entries.append(entry)
        return key

    def reflection_for_root(self, root):
        beta = np.asarray(root, dtype=np.float64).reshape(3)
        norm = float(beta.T @ self.gram @ beta)
        if abs(norm) <= self.tolerance:
            raise ValueError("An isotropic vector cannot define a Coxeter reflection.")
        return np.eye(3) - (2.0 / norm) * np.outer(beta, beta.T @ self.gram)

    def _record_root(self, root, depth):
        positive = self._positive_root(root)
        key = self.root_key(positive)
        record = self._root_registry.get(key)
        if record is None:
            record = RootRecord(
                key=key,
                coefficients=positive,
                reflection=self.reflection_for_root(positive),
                depth=int(depth),
            )
            self._root_registry[key] = record
        else:
            record.depth = min(record.depth, int(depth))
        return record

    def enumerate_chambers(self, max_length):
        if not isinstance(max_length, int) or max_length < 0:
            raise ValueError("max_length must be a nonnegative integer.")

        identity = np.eye(3)
        identity_key = self.element_key(identity)
        chambers = {
            identity_key: ChamberRecord(
                key=identity_key,
                word="e",
                matrix=identity.copy(),
                aliases={"e"},
                length=0,
            )
        }
        queue = deque([identity_key])

        while queue:
            chamber_key = queue.popleft()
            chamber = chambers[chamber_key]
            word_digits = "" if chamber.word == "e" else chamber.word
            matrix = chamber.matrix
            matrix_inverse = np.linalg.inv(matrix)

            for index, generator in enumerate(GENERATOR_ORDER):
                candidate_word = word_digits + generator
                candidate_matrix = matrix @ self.generator_matrices[generator]
                candidate_key = self.element_key(candidate_matrix)
                wall_root = matrix @ self.simple_roots[index]
                root_record = self._record_root(
                    wall_root,
                    depth=chamber.length + 1,
                )
                wall_reflection = (
                    matrix
                    @ self.generator_matrices[generator]
                    @ matrix_inverse
                )
                chamber.walls[generator] = WallRecord(
                    root_key=root_record.key,
                    reflection=wall_reflection,
                    wall_type=generator,
                    chamber_key=chamber.key,
                    neighbor_key=candidate_key,
                )

                existing = chambers.get(candidate_key)
                if existing is not None:
                    existing.aliases.add(candidate_word or "e")
                    continue
                if chamber.length >= max_length:
                    continue

                candidate = ChamberRecord(
                    key=candidate_key,
                    word=candidate_word,
                    matrix=candidate_matrix.copy(),
                    aliases={candidate_word},
                    length=chamber.length + 1,
                    inversions=chamber.inversions + (root_record.key,),
                )
                chambers[candidate_key] = candidate
                queue.append(candidate_key)

        return chambers

    def canonicalize_word(self, word):
        normalized = self._normalize_word(word)
        display_word = normalized or "e"
        chambers = self.enumerate_chambers(len(normalized))
        key = self.element_key(self.word_matrix(normalized))
        chamber = chambers.get(key)
        if chamber is None:
            target = self.word_matrix(normalized)
            chamber = next(
                (
                    record
                    for record in chambers.values()
                    if np.allclose(
                        self.word_matrix(record.word),
                        target,
                        atol=self.tolerance * 10,
                        rtol=self.tolerance * 10,
                    )
                ),
                None,
            )
        if chamber is None:
            raise ArithmeticError(f"Cannot find the chamber corresponding to the word {word}.")
        chamber.aliases.add(display_word)
        return chamber

    def enumerate_small_roots(self):
        small = {}
        queue = deque()

        for simple_root in self.simple_roots:
            record = self._record_root(simple_root, depth=1)
            if record.key not in small:
                record.is_small = True
                small[record.key] = record
                queue.append(record)

        while queue:
            record = queue.popleft()
            beta = record.coefficients
            for index, generator in enumerate(GENERATOR_ORDER):
                inner_product = float(
                    beta.T @ self.gram @ self.simple_roots[index]
                )
                if not (-1.0 + self.tolerance < inner_product < -self.tolerance):
                    continue
                reflected = self.generator_matrices[generator] @ beta
                candidate = self._record_root(reflected, depth=record.depth + 1)
                if candidate.key in small:
                    continue
                candidate.is_small = True
                small[candidate.key] = candidate
                queue.append(candidate)

        return small

    def enumerate_m_small_roots(self, m):
        """Enumerate Sigma_m using the dominance-depth reflection rule."""
        if not isinstance(m, int) or m < 0:
            raise ValueError("m must be a nonnegative integer.")
        if m == 0:
            return self.enumerate_small_roots()

        roots = dict(self.enumerate_small_roots())
        dominance_depths = {key: 0 for key in roots}
        queue = deque(roots)

        while queue:
            key = queue.popleft()
            record = roots[key]
            beta = record.coefficients
            current_depth = dominance_depths[key]
            for index, generator in enumerate(GENERATOR_ORDER):
                reflected = self.generator_matrices[generator] @ beta
                if np.all(reflected <= self.tolerance * 10):
                    continue
                reflected = self._positive_root(reflected)
                candidate_key = self.root_key(reflected)
                inner_product = float(
                    beta.T @ self.gram @ self.simple_roots[index]
                )
                comparison_tolerance = self.tolerance * 100
                if inner_product >= 1.0 - comparison_tolerance:
                    candidate_depth = current_depth - 1
                elif inner_product <= -1.0 + comparison_tolerance:
                    candidate_depth = current_depth + 1
                else:
                    candidate_depth = current_depth

                if candidate_depth < 0 or candidate_depth > m:
                    continue
                known_depth = dominance_depths.get(candidate_key)
                if known_depth is not None:
                    if known_depth != candidate_depth:
                        raise ArithmeticError(
                            "Inconsistent dominance depths were obtained for the same root."
                        )
                    continue

                candidate = self._record_root(
                    reflected,
                    depth=record.depth + 1,
                )
                roots[candidate_key] = candidate
                dominance_depths[candidate_key] = candidate_depth
                queue.append(candidate_key)

        return roots

    def cone_contains(self, root, generators):
        target = np.asarray(root, dtype=np.float64).reshape(3)
        candidates = [
            np.asarray(generator, dtype=np.float64).reshape(3)
            for generator in generators
        ]
        if np.linalg.norm(target) <= self.tolerance:
            return True
        if not candidates:
            return False

        residual_tolerance = max(self.tolerance * 100, 1e-8) * (
            1.0 + np.linalg.norm(target)
        )
        for size in range(1, min(3, len(candidates)) + 1):
            for selected in combinations(candidates, size):
                matrix = np.column_stack(selected)
                coefficients, _, _, _ = np.linalg.lstsq(matrix, target, rcond=None)
                residual = np.linalg.norm(matrix @ coefficients - target)
                if (
                    residual <= residual_tolerance
                    and np.all(coefficients >= -self.tolerance * 100)
                ):
                    return True
        return False

    def is_low_element(self, chamber, small_roots=None):
        small = small_roots or self.enumerate_small_roots()
        if not chamber.inversions:
            return True

        inversion_records = [
            self._root_registry[root_key]
            for root_key in chamber.inversions
        ]
        small_inversions = [
            record.coefficients
            for record in inversion_records
            if record.key in small
        ]
        return all(
            self.cone_contains(record.coefficients, small_inversions)
            for record in inversion_records
        )

    def enumerate_low_elements(self):
        """Enumerate L_0 using its suffix-closed Garside-shadow property.

        Every non-identity low element has a reduced expression whose suffix
        after deleting the first generator is again low.  Starting at the
        identity and prepending generators therefore terminates exactly when
        the finite set L_0 is exhausted.
        """

        small_roots = self.enumerate_small_roots()
        identity = self.canonicalize_word("e")
        low_elements = {identity.key: identity}
        frontier = [identity]

        while frontier:
            next_frontier = []
            for chamber in frontier:
                suffix = "" if chamber.word == "e" else chamber.word
                for generator in GENERATOR_ORDER:
                    candidate = self.canonicalize_word(generator + suffix)
                    if candidate.length != chamber.length + 1:
                        continue
                    if candidate.key in low_elements:
                        continue
                    if not self.is_low_element(candidate, small_roots):
                        continue
                    low_elements[candidate.key] = candidate
                    next_frontier.append(candidate)
            frontier = next_frontier

        return low_elements

    def enumerate_m_low_elements(self, m):
        """Enumerate the finite Garside shadow L_m."""
        if not isinstance(m, int) or m < 0:
            raise ValueError("m must be a nonnegative integer.")
        if m == 0:
            return self.enumerate_low_elements()

        m_small_roots = self.enumerate_m_small_roots(m)
        identity = self.canonicalize_word("e")
        low_elements = {identity.key: identity}
        frontier = [identity]

        while frontier:
            next_frontier = []
            for chamber in frontier:
                suffix = "" if chamber.word == "e" else chamber.word
                for generator in GENERATOR_ORDER:
                    candidate = self.canonicalize_word(generator + suffix)
                    if candidate.length != chamber.length + 1:
                        continue
                    if candidate.key in low_elements:
                        continue
                    if not self.is_low_element(candidate, m_small_roots):
                        continue
                    low_elements[candidate.key] = candidate
                    next_frontier.append(candidate)
            frontier = next_frontier

        return low_elements
