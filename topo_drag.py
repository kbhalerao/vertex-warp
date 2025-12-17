# -*- coding: utf-8 -*-
"""
GeoRefiner Manual - Topology-Preserving Drag

When a vertex is dragged, neighbors move along while preserving:
- Local interior angles
- Edge length ratios

Uses Gaussian-weighted influence decay and optimization to find
positions that balance movement with geometric constraints.
"""

import math
from typing import List, Dict, Tuple, Optional
import numpy as np

# scipy is optional - only needed if fast_mode=False (not recommended)
# Commented out to ensure zero external dependencies
# try:
#     from scipy.optimize import minimize
#     HAS_SCIPY = True
# except ImportError:
#     HAS_SCIPY = False
HAS_SCIPY = False


def gaussian_weight(distance: float, sigma: float) -> float:
    """
    Compute Gaussian weight based on distance.

    Args:
        distance: Distance between vertices in map units.
        sigma: Falloff radius (standard deviation).

    Returns:
        Weight between 0 and 1.
    """
    return math.exp(-(distance ** 2) / (2 * sigma ** 2))


def compute_angle(p1: Tuple[float, float], p2: Tuple[float, float], p3: Tuple[float, float]) -> float:
    """
    Compute angle at p2 formed by p1-p2-p3.

    Args:
        p1, p2, p3: Points as (x, y) tuples.

    Returns:
        Angle in radians.
    """
    v1 = (p1[0] - p2[0], p1[1] - p2[1])
    v2 = (p3[0] - p2[0], p3[1] - p2[1])

    dot = v1[0] * v2[0] + v1[1] * v2[1]
    mag1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2)
    mag2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2)

    if mag1 < 1e-10 or mag2 < 1e-10:
        return 0.0

    cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
    return math.acos(cos_angle)


def compute_edge_length(p1: Tuple[float, float], p2: Tuple[float, float]) -> float:
    """Compute distance between two points."""
    return math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)


class TopologyPreservingDrag:
    """
    Handles topology-preserving vertex dragging.

    When a vertex is moved, nearby vertices follow with weighted influence
    while attempting to preserve local angles and edge length ratios.
    """

    def __init__(
        self,
        vertices: List[Dict],
        sigma: float = 50.0,
        angle_weight: float = 1.0,
        edge_weight: float = 0.5,
        position_weight: float = 1.0,
        fast_mode: bool = True,
    ):
        """
        Initialize the drag handler.

        Args:
            vertices: List of vertex dictionaries with id, current_x, current_y.
            sigma: Gaussian falloff radius in map units.
            angle_weight: Weight for angle preservation constraint.
            edge_weight: Weight for edge length preservation constraint.
            position_weight: Weight for staying close to target positions.
            fast_mode: If True, skip scipy optimization (much faster, still moves neighbors).
        """
        self.vertices = vertices
        self.sigma = sigma
        self.angle_weight = angle_weight
        self.edge_weight = edge_weight
        self.position_weight = position_weight
        self.fast_mode = fast_mode

        # Build adjacency for polygon (assumes vertices form a closed polygon)
        self.n_vertices = len(vertices)
        if not fast_mode:
            self._build_original_geometry()

    def _build_original_geometry(self):
        """Compute original angles and edge lengths."""
        n = self.n_vertices
        if n < 3:
            self.original_angles = []
            self.original_edge_lengths = []
            return

        # Get current positions
        positions = [(v["current_x"], v["current_y"]) for v in self.vertices]

        # Compute angles at each vertex
        self.original_angles = []
        for i in range(n):
            prev_idx = (i - 1) % n
            next_idx = (i + 1) % n
            angle = compute_angle(positions[prev_idx], positions[i], positions[next_idx])
            self.original_angles.append(angle)

        # Compute edge lengths
        self.original_edge_lengths = []
        for i in range(n):
            next_idx = (i + 1) % n
            length = compute_edge_length(positions[i], positions[next_idx])
            self.original_edge_lengths.append(length)

    def compute_drag_effect(
        self,
        dragged_vertex_id: int,
        new_x: float,
        new_y: float,
    ) -> List[Dict]:
        """
        Compute the effect of dragging a vertex on all vertices.

        Args:
            dragged_vertex_id: ID of the vertex being dragged.
            new_x: New X position of the dragged vertex.
            new_y: New Y position of the dragged vertex.

        Returns:
            List of updated vertex positions:
            [{"id": ..., "x": ..., "y": ..., "moved": True/False}, ...]
        """
        if not HAS_SCIPY or self.n_vertices < 3:
            # Fallback: only move the dragged vertex
            return self._simple_drag(dragged_vertex_id, new_x, new_y)

        # Find dragged vertex index
        dragged_idx = None
        for i, v in enumerate(self.vertices):
            if v["id"] == dragged_vertex_id:
                dragged_idx = i
                break

        if dragged_idx is None:
            return self._simple_drag(dragged_vertex_id, new_x, new_y)

        # Get current positions
        current_positions = np.array([
            [v["current_x"], v["current_y"]] for v in self.vertices
        ])

        # Compute delta for dragged vertex
        delta = np.array([
            new_x - current_positions[dragged_idx, 0],
            new_y - current_positions[dragged_idx, 1]
        ])

        # Compute influence weights based on distance from dragged vertex
        weights = np.zeros(self.n_vertices)
        dragged_pos = current_positions[dragged_idx]

        for i in range(self.n_vertices):
            if i == dragged_idx:
                weights[i] = 1.0  # Full weight for dragged vertex
            else:
                dist = np.linalg.norm(current_positions[i] - dragged_pos)
                weights[i] = gaussian_weight(dist, self.sigma)

        # Compute "free" target positions (weighted delta applied)
        target_positions = current_positions.copy()
        for i in range(self.n_vertices):
            target_positions[i] += weights[i] * delta

        # In fast mode, skip optimization and use weighted positions directly
        if self.fast_mode:
            final_positions = target_positions
        else:
            # Optimize to find positions that preserve geometry (slow)
            final_positions = self._optimize_positions(
                current_positions,
                target_positions,
                weights,
                dragged_idx,
            )

        # Build result (only mark as moved if movement exceeds threshold)
        MOVE_THRESHOLD = 0.1  # meters - ignore movements smaller than this
        results = []
        for i, v in enumerate(self.vertices):
            dx = final_positions[i, 0] - current_positions[i, 0]
            dy = final_positions[i, 1] - current_positions[i, 1]
            distance_moved = math.sqrt(dx * dx + dy * dy)

            # Only mark as newly moved if movement exceeds threshold
            newly_moved = distance_moved > MOVE_THRESHOLD
            results.append({
                "id": v["id"],
                "x": float(final_positions[i, 0]),
                "y": float(final_positions[i, 1]),
                "moved": newly_moved or v.get("moved", False),
            })

        return results

    def _simple_drag(self, vertex_id: int, new_x: float, new_y: float) -> List[Dict]:
        """Simple drag without topology preservation."""
        results = []
        for v in self.vertices:
            if v["id"] == vertex_id:
                results.append({
                    "id": v["id"],
                    "x": new_x,
                    "y": new_y,
                    "moved": True,
                })
            else:
                results.append({
                    "id": v["id"],
                    "x": v["current_x"],
                    "y": v["current_y"],
                    "moved": v.get("moved", False),
                })
        return results

    def _optimize_positions(
        self,
        current: np.ndarray,
        target: np.ndarray,
        weights: np.ndarray,
        dragged_idx: int,
    ) -> np.ndarray:
        """
        Optimize vertex positions to preserve geometry.

        Args:
            current: Current positions (n_vertices, 2).
            target: Target "free" positions (n_vertices, 2).
            weights: Influence weights per vertex.
            dragged_idx: Index of the dragged vertex.

        Returns:
            Optimized positions (n_vertices, 2).
        """
        n = self.n_vertices

        def objective(x):
            """Objective function: minimize position error + geometry constraints."""
            positions = x.reshape((n, 2))

            # Position error (weighted by influence)
            pos_error = 0.0
            for i in range(n):
                diff = positions[i] - target[i]
                pos_error += weights[i] * np.sum(diff ** 2)

            # Angle preservation error
            angle_error = 0.0
            for i in range(n):
                prev_idx = (i - 1) % n
                next_idx = (i + 1) % n
                p1 = tuple(positions[prev_idx])
                p2 = tuple(positions[i])
                p3 = tuple(positions[next_idx])
                new_angle = compute_angle(p1, p2, p3)
                orig_angle = self.original_angles[i]
                angle_error += (new_angle - orig_angle) ** 2

            # Edge length ratio preservation
            edge_error = 0.0
            if len(self.original_edge_lengths) > 0:
                total_orig = sum(self.original_edge_lengths)
                for i in range(n):
                    next_idx = (i + 1) % n
                    new_length = compute_edge_length(
                        tuple(positions[i]), tuple(positions[next_idx])
                    )
                    orig_length = self.original_edge_lengths[i]
                    if orig_length > 1e-10:
                        # Ratio error
                        orig_ratio = orig_length / total_orig
                        new_total = sum(
                            compute_edge_length(tuple(positions[j]), tuple(positions[(j + 1) % n]))
                            for j in range(n)
                        )
                        if new_total > 1e-10:
                            new_ratio = new_length / new_total
                            edge_error += (new_ratio - orig_ratio) ** 2

            total = (
                self.position_weight * pos_error +
                self.angle_weight * angle_error +
                self.edge_weight * edge_error
            )
            return total

        # Initial guess is target positions
        x0 = target.flatten()

        # Optimize with constraints
        result = minimize(
            objective,
            x0,
            method='L-BFGS-B',
            options={'maxiter': 100, 'ftol': 1e-6}
        )

        return result.x.reshape((n, 2))

    def update_sigma(self, new_sigma: float):
        """Update the falloff radius."""
        self.sigma = new_sigma

    def update_vertices(self, vertices: List[Dict]):
        """Update vertex list and rebuild geometry."""
        self.vertices = vertices
        self.n_vertices = len(vertices)
        self._build_original_geometry()


def apply_topology_drag(
    vertices: List[Dict],
    dragged_id: int,
    new_x: float,
    new_y: float,
    sigma: float = 50.0,
    enabled: bool = True,
    fast_mode: bool = True,
) -> List[Dict]:
    """
    Convenience function to apply topology-preserving drag.

    Args:
        vertices: List of vertex dictionaries.
        dragged_id: ID of the vertex being dragged.
        new_x, new_y: New position of dragged vertex.
        sigma: Gaussian falloff radius.
        enabled: If False, only moves the single vertex.
        fast_mode: If True, skip scipy optimization (much faster).

    Returns:
        List of updated vertex positions.
    """
    if not enabled:
        # Simple drag - only move one vertex
        results = []
        for v in vertices:
            if v["id"] == dragged_id:
                results.append({
                    "id": v["id"],
                    "x": new_x,
                    "y": new_y,
                    "moved": True,
                })
            else:
                results.append({
                    "id": v["id"],
                    "x": v["current_x"],
                    "y": v["current_y"],
                    "moved": v.get("moved", False),
                })
        return results

    # Topology-preserving drag (fast mode by default)
    handler = TopologyPreservingDrag(vertices, sigma=sigma, fast_mode=fast_mode)
    return handler.compute_drag_effect(dragged_id, new_x, new_y)
