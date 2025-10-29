"""Portfolio optimisation analysis for MATH2070 assignment tasks (excluding Question 6).

The script implements all calculations using the Python standard library only so that it can
run in the execution environment without external dependencies such as NumPy or SciPy.
It reproduces the requested figures and numerical answers for Questions 1-5 and 7.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

DATA_FILE = Path("project_data.csv")
FIGURE_DIR = Path("figures")


@dataclass
class MarketData:
    dates: List[str]
    assets: List[str]
    prices: Dict[str, List[float]]
    returns: Dict[str, List[float]]
    mean_vector: List[float]
    covariance: List[List[float]]


def read_market_data(path: Path) -> MarketData:
    """Load price data and compute daily simple returns, means and covariance matrix."""

    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    dates = [row[next(iter(row))] for row in rows]
    assets = [name for name in rows[0].keys() if name not in {"Date", "\ufeffDate"}]

    prices: Dict[str, List[float]] = {asset: [] for asset in assets}
    for row in rows:
        for asset in assets:
            prices[asset].append(float(row[asset]))

    returns: Dict[str, List[float]] = {asset: [] for asset in assets}
    for asset in assets:
        series = prices[asset]
        returns[asset] = [series[i] / series[i - 1] - 1.0 for i in range(1, len(series))]

    count = len(next(iter(returns.values())))
    mean_vector = [sum(returns[asset]) / count for asset in assets]

    covariance = [[0.0 for _ in assets] for _ in assets]
    for i, asset_i in enumerate(assets):
        mean_i = mean_vector[i]
        for j, asset_j in enumerate(assets):
            mean_j = mean_vector[j]
            covariance[i][j] = (
                sum(
                    (returns[asset_i][k] - mean_i)
                    * (returns[asset_j][k] - mean_j)
                    for k in range(count)
                )
                / (count - 1)
            )

    return MarketData(dates, assets, prices, returns, mean_vector, covariance)


def gaussian_elimination_solve(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> List[float]:
    """Solve a linear system using Gaussian elimination with partial pivoting."""

    n = len(matrix)
    aug = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]

    for col in range(n):
        pivot_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
        pivot_val = aug[pivot_row][col]
        if abs(pivot_val) < 1e-15:
            raise ValueError("Singular matrix encountered during elimination")
        if pivot_row != col:
            aug[col], aug[pivot_row] = aug[pivot_row], aug[col]

        pivot_val = aug[col][col]
        for j in range(col, n + 1):
            aug[col][j] /= pivot_val

        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            if abs(factor) < 1e-15:
                continue
            for j in range(col, n + 1):
                aug[row][j] -= factor * aug[col][j]

    return [aug[i][n] for i in range(n)]


def dot(a: Sequence[float], b: Sequence[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def mat_vec_mul(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> List[float]:
    return [dot(row, vector) for row in matrix]


def portfolio_statistics(weights: Sequence[float], mean_vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> Tuple[float, float, float]:
    mean_return = dot(weights, mean_vector)
    variance = dot(weights, mat_vec_mul(covariance, weights))
    return mean_return, variance, math.sqrt(max(variance, 0.0))


def constrained_optimal_weights(
    covariance: Sequence[Sequence[float]],
    mean_vector: Sequence[float],
    risk_tolerance: float,
) -> Tuple[List[float], float, float]:
    """Unconstrained Markowitz optimiser with only the budget constraint."""

    ones = [1.0] * len(mean_vector)
    c_inv_r = gaussian_elimination_solve(covariance, mean_vector)
    c_inv_ones = gaussian_elimination_solve(covariance, ones)
    alpha = dot(ones, c_inv_r)
    beta = dot(ones, c_inv_ones)

    weights = [
        risk_tolerance * (c_inv_r[i] - (alpha / beta) * c_inv_ones[i])
        + c_inv_ones[i] / beta
        for i in range(len(mean_vector))
    ]
    mean_return, variance, std_dev = portfolio_statistics(weights, mean_vector, covariance)
    return weights, mean_return, std_dev


def global_minimum_variance_weights(
    covariance: Sequence[Sequence[float]],
) -> List[float]:
    ones = [1.0] * len(covariance)
    c_inv_ones = gaussian_elimination_solve(covariance, ones)
    denom = dot(ones, c_inv_ones)
    return [val / denom for val in c_inv_ones]


def efficient_frontier(
    covariance: Sequence[Sequence[float]],
    mean_vector: Sequence[float],
    steps: int = 100,
    mu_scale: float = 1.6,
) -> List[Tuple[float, float]]:
    ones = [1.0] * len(mean_vector)
    c_inv_r = gaussian_elimination_solve(covariance, mean_vector)
    c_inv_ones = gaussian_elimination_solve(covariance, ones)
    A = dot(ones, c_inv_ones)
    B = dot(ones, c_inv_r)
    C = dot(mean_vector, c_inv_r)
    D = A * C - B * B

    mu_min = dot(mean_vector, global_minimum_variance_weights(covariance))
    mu_max = max(mean_vector) * mu_scale
    mu_values = [mu_min + (mu_max - mu_min) * i / (steps - 1) for i in range(steps)]

    frontier = []
    for mu in mu_values:
        variance = (A * mu * mu - 2 * B * mu + C) / D
        frontier.append((math.sqrt(max(variance, 0.0)), mu))
    return frontier


def enumerate_no_short_solutions(
    covariance: Sequence[Sequence[float]],
    mean_vector: Sequence[float],
    risk_tolerance: float,
) -> Tuple[List[float], float, float]:
    n = len(mean_vector)
    best_weights: List[float] | None = None
    best_objective = float("-inf")
    best_stats: Tuple[float, float, float] = (0.0, 0.0, 0.0)

    indices = range(n)
    for r in range(1, n + 1):
        for combo in combinations(indices, r):
            sub_cov = [[covariance[i][j] for j in combo] for i in combo]
            sub_mean = [mean_vector[i] for i in combo]
            ones = [1.0] * len(combo)
            try:
                c_inv_r = gaussian_elimination_solve(sub_cov, sub_mean)
                c_inv_ones = gaussian_elimination_solve(sub_cov, ones)
            except ValueError:
                continue

            alpha = dot(ones, c_inv_r)
            beta = dot(ones, c_inv_ones)
            weights_subset = [
                risk_tolerance * (c_inv_r[i] - (alpha / beta) * c_inv_ones[i])
                + c_inv_ones[i] / beta
                for i in range(len(combo))
            ]

            if any(weight < -1e-9 for weight in weights_subset):
                continue

            weights = [0.0] * n
            for idx, weight in zip(combo, weights_subset):
                weights[idx] = weight

            if abs(sum(weights) - 1.0) > 1e-6:
                continue

            mean_return, variance, std_dev = portfolio_statistics(weights, mean_vector, covariance)
            objective = risk_tolerance * mean_return - 0.5 * variance
            if objective > best_objective:
                best_objective = objective
                best_weights = weights
                best_stats = (mean_return, variance, std_dev)

    if best_weights is None:
        raise RuntimeError("No feasible portfolio found without short selling")

    return best_weights, best_stats[0], best_stats[2]


def kkt_residuals(
    covariance: Sequence[Sequence[float]],
    mean_vector: Sequence[float],
    weights: Sequence[float],
    risk_tolerance: float,
) -> Tuple[float, List[float]]:
    covariance_times_weights = mat_vec_mul(covariance, weights)
    lambdas = [risk_tolerance * mean_vector[i] - covariance_times_weights[i] for i, w in enumerate(weights) if w > 1e-8]
    lam = sum(lambdas) / len(lambdas)
    gamma = [covariance_times_weights[i] - risk_tolerance * mean_vector[i] + lam for i in range(len(weights))]
    return lam, gamma


def tangent_portfolio(
    covariance: Sequence[Sequence[float]],
    mean_vector: Sequence[float],
    risk_free_rate: float,
) -> Tuple[List[float], float, float]:
    excess = [mu - risk_free_rate for mu in mean_vector]
    c_inv_excess = gaussian_elimination_solve(covariance, excess)
    denom = sum(c_inv_excess)
    weights = [val / denom for val in c_inv_excess]
    mean_return, variance, std_dev = portfolio_statistics(weights, mean_vector, covariance)
    return weights, mean_return, std_dev


def risk_free_optimal_portfolio(
    covariance: Sequence[Sequence[float]],
    mean_vector: Sequence[float],
    risk_tolerance: float,
    risk_free_rate: float,
) -> Tuple[List[float], float, float, float]:
    excess = [mu - risk_free_rate for mu in mean_vector]
    c_inv_excess = gaussian_elimination_solve(covariance, excess)
    risky_weights = [risk_tolerance * val for val in c_inv_excess]
    risky_mean = dot(risky_weights, excess)
    risky_variance = dot(risky_weights, mat_vec_mul(covariance, risky_weights))
    portfolio_mean = risk_free_rate + risky_mean
    portfolio_std = math.sqrt(max(risky_variance, 0.0))
    risk_free_weight = 1.0 - sum(risky_weights)
    return risky_weights, risk_free_weight, portfolio_mean, portfolio_std


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def build_svg_axes(
    width: int,
    height: int,
    margin: int,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
) -> List[str]:
    elements: List[str] = []
    # axes lines
    elements.append(
        f'<line x1="{margin}" y1="{height - margin}" x2="{width - margin}" y2="{height - margin}" '
        'stroke="black" stroke-width="2" />'
    )
    elements.append(
        f'<line x1="{margin}" y1="{height - margin}" x2="{margin}" y2="{margin}" stroke="black" stroke-width="2" />'
    )

    # x ticks
    for frac in range(0, 11):
        sigma = x_min + (x_max - x_min) * frac / 10
        x = margin + (sigma - x_min) / (x_max - x_min) * (width - 2 * margin)
        elements.append(
            f'<line x1="{x}" y1="{height - margin}" x2="{x}" y2="{height - margin + 8}" stroke="black" stroke-width="1" />'
        )
        elements.append(
            f'<text x="{x}" y="{height - margin + 24}" font-size="14" text-anchor="middle">{sigma:.3f}</text>'
        )

    # y ticks
    for frac in range(0, 11):
        mu = y_min + (y_max - y_min) * frac / 10
        y = height - margin - (mu - y_min) / (y_max - y_min) * (height - 2 * margin)
        elements.append(
            f'<line x1="{margin}" y1="{y}" x2="{margin - 8}" y2="{y}" stroke="black" stroke-width="1" />'
        )
        elements.append(
            f'<text x="{margin - 12}" y="{y + 5}" font-size="14" text-anchor="end">{mu:.4f}</text>'
        )

    elements.append(
        f'<text x="{(width) / 2}" y="{height - margin + 50}" font-size="16" text-anchor="middle">Portfolio risk σ</text>'
    )
    elements.append(
        f'<text x="{margin - 60}" y="{(height) / 2}" font-size="16" text-anchor="middle" transform="rotate(-90 {margin - 60},{height / 2})">'
        'Expected return μ</text>'
    )

    return elements


def point_to_svg(
    sigma: float,
    mu: float,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    width: int,
    height: int,
    margin: int,
) -> Tuple[float, float]:
    x = margin + (sigma - x_min) / (x_max - x_min) * (width - 2 * margin)
    y = height - margin - (mu - y_min) / (y_max - y_min) * (height - 2 * margin)
    return x, y


def add_polyline(points: Sequence[Tuple[float, float]], **attrs: str) -> str:
    attr_str = " ".join(f"{key}='{value}'" for key, value in attrs.items())
    coords = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return f"<polyline points='{coords}' {attr_str} />"


def generate_question4_svg(
    data: MarketData,
    risk_tolerance: float,
    weights: Sequence[float],
    gmvp_weights: Sequence[float],
    frontier: Sequence[Tuple[float, float]],
    file_path: Path,
) -> None:
    width, height, margin = 900, 650, 90
    x_min, x_max = 0.005, 0.02
    y_min, y_max = 0.0004, 0.0015

    svg_elements = [
        "<?xml version='1.0' encoding='UTF-8'?>",
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
    ]
    svg_elements.extend(build_svg_axes(width, height, margin, x_min, x_max, y_min, y_max))

    # Minimum variance frontier (entire)
    frontier_points = [
        point_to_svg(sigma, mu, x_min, x_max, y_min, y_max, width, height, margin)
        for sigma, mu in frontier
    ]
    svg_elements.append(
        add_polyline(frontier_points, fill="none", stroke="#6baed6", **{"stroke-width": "2", "stroke-dasharray": "6 4"})
    )

    # Efficient frontier (upper half)
    mu_gmv = portfolio_statistics(gmvp_weights, data.mean_vector, data.covariance)[0]
    efficient_points = [
        pt
        for (sigma, mu), pt in zip(frontier, frontier_points)
        if mu >= mu_gmv - 1e-9
    ]
    svg_elements.append(
        add_polyline(efficient_points, fill="none", stroke="#08519c", **{"stroke-width": "3"})
    )

    # Plot assets
    for idx, asset in enumerate(data.assets):
        sigma = math.sqrt(data.covariance[idx][idx])
        mu = data.mean_vector[idx]
        x, y = point_to_svg(sigma, mu, x_min, x_max, y_min, y_max, width, height, margin)
        svg_elements.append(
            f"<circle cx='{x:.2f}' cy='{y:.2f}' r='5' fill='#de2d26' />"
        )
        svg_elements.append(
            f"<text x='{x + 8:.2f}' y='{y - 8:.2f}' font-size='14' fill='#de2d26'>{asset}</text>"
        )

    # Optimal unrestricted portfolio
    mu_opt, _, sigma_opt = portfolio_statistics(weights, data.mean_vector, data.covariance)
    x_opt, y_opt = point_to_svg(sigma_opt, mu_opt, x_min, x_max, y_min, y_max, width, height, margin)
    svg_elements.append(
        f"<rect x='{x_opt - 6:.2f}' y='{y_opt - 6:.2f}' width='12' height='12' fill='#238b45' />"
    )
    svg_elements.append(
        f"<text x='{x_opt + 10:.2f}' y='{y_opt + 4:.2f}' font-size='14' fill='#238b45'>Optimal Q3</text>"
    )

    # GMVP
    mu_gmv, _, sigma_gmv = portfolio_statistics(gmvp_weights, data.mean_vector, data.covariance)
    x_gmv, y_gmv = point_to_svg(sigma_gmv, mu_gmv, x_min, x_max, y_min, y_max, width, height, margin)
    svg_elements.append(
        f"<polygon points='{x_gmv - 6:.2f},{y_gmv:.2f} {x_gmv:.2f},{y_gmv - 8:.2f} {x_gmv + 6:.2f},{y_gmv:.2f}' fill='#756bb1' />"
    )
    svg_elements.append(
        f"<text x='{x_gmv + 10:.2f}' y='{y_gmv + 4:.2f}' font-size='14' fill='#756bb1'>GMV</text>"
    )

    svg_elements.append("</svg>")
    file_path.write_text("\n".join(svg_elements))


def generate_question7_svg(
    data: MarketData,
    frontier: Sequence[Tuple[float, float]],
    q3_weights: Sequence[float],
    q7_weights: Sequence[float],
    q7_risk_free_weight: float,
    tangent_weights: Sequence[float],
    risk_free_rate: float,
    market_point: Tuple[float, float, float],
    file_path: Path,
) -> None:
    width, height, margin = 900, 650, 90
    x_min, x_max = 0.0, 0.02
    y_min, y_max = 0.0, 0.0015

    svg_elements = [
        "<?xml version='1.0' encoding='UTF-8'?>",
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
    ]
    svg_elements.extend(build_svg_axes(width, height, margin, x_min, x_max, y_min, y_max))

    # Original frontier
    frontier_points = [
        point_to_svg(sigma, mu, x_min, x_max, y_min, y_max, width, height, margin)
        for sigma, mu in frontier
    ]
    svg_elements.append(
        add_polyline(frontier_points, fill="none", stroke="#6baed6", **{"stroke-width": "2", "stroke-dasharray": "6 4"})
    )

    # Efficient half
    mu_gmv = portfolio_statistics(global_minimum_variance_weights(data.covariance), data.mean_vector, data.covariance)[0]
    efficient_points = [
        pt
        for (sigma, mu), pt in zip(frontier, frontier_points)
        if mu >= mu_gmv - 1e-9
    ]
    svg_elements.append(
        add_polyline(efficient_points, fill="none", stroke="#08519c", **{"stroke-width": "3"})
    )

    # Assets
    for idx, asset in enumerate(data.assets):
        sigma = math.sqrt(data.covariance[idx][idx])
        mu = data.mean_vector[idx]
        x, y = point_to_svg(sigma, mu, x_min, x_max, y_min, y_max, width, height, margin)
        svg_elements.append(f"<circle cx='{x:.2f}' cy='{y:.2f}' r='5' fill='#de2d26' />")
        svg_elements.append(
            f"<text x='{x + 8:.2f}' y='{y - 8:.2f}' font-size='14' fill='#de2d26'>{asset}</text>"
        )

    # Risk-free point
    x_rf, y_rf = point_to_svg(0.0, risk_free_rate, x_min, x_max, y_min, y_max, width, height, margin)
    svg_elements.append(f"<circle cx='{x_rf:.2f}' cy='{y_rf:.2f}' r='5' fill='#feb24c' />")
    svg_elements.append(
        f"<text x='{x_rf + 10:.2f}' y='{y_rf + 4:.2f}' font-size='14' fill='#feb24c'>Risk-free</text>"
    )

    # Capital market line
    slope = (market_point[1] - risk_free_rate) / market_point[2]
    line_points = []
    for sigma in [0.0, 0.02]:
        mu = risk_free_rate + slope * sigma
        line_points.append(
            point_to_svg(sigma, mu, x_min, x_max, y_min, y_max, width, height, margin)
        )
    svg_elements.append(add_polyline(line_points, stroke="#31a354", fill="none", **{"stroke-width": "2"}))
    svg_elements.append(
        f"<text x='{(line_points[0][0] + line_points[1][0]) / 2:.2f}' y='{(line_points[0][1] + line_points[1][1]) / 2 - 12:.2f}' font-size='14' fill='#31a354'>CML</text>"
    )

    # Original optimal portfolio
    mu_q3, _, sigma_q3 = portfolio_statistics(q3_weights, data.mean_vector, data.covariance)
    x_q3, y_q3 = point_to_svg(sigma_q3, mu_q3, x_min, x_max, y_min, y_max, width, height, margin)
    svg_elements.append(f"<rect x='{x_q3 - 6:.2f}' y='{y_q3 - 6:.2f}' width='12' height='12' fill='#238b45' />")
    svg_elements.append(
        f"<text x='{x_q3 + 10:.2f}' y='{y_q3 + 4:.2f}' font-size='14' fill='#238b45'>Optimal Q3</text>"
    )

    # Optimal portfolio with risk-free asset
    mu_q7, _, sigma_q7 = portfolio_statistics(q7_weights, data.mean_vector, data.covariance)
    mu_q7_total = risk_free_rate + dot(q7_weights, [mu - risk_free_rate for mu in data.mean_vector])
    x_q7, y_q7 = point_to_svg(sigma_q7, mu_q7_total, x_min, x_max, y_min, y_max, width, height, margin)
    svg_elements.append(f"<polygon points='{x_q7 - 6:.2f},{y_q7 + 6:.2f} {x_q7 + 6:.2f},{y_q7 + 6:.2f} {x_q7:.2f},{y_q7 - 6:.2f}' fill='#756bb1' />")
    svg_elements.append(
        f"<text x='{x_q7 + 10:.2f}' y='{y_q7 + 4:.2f}' font-size='14' fill='#756bb1'>Optimal Q7a</text>"
    )

    # Market / tangency portfolio
    mu_market, _, sigma_market = market_point
    x_market, y_market = point_to_svg(sigma_market, mu_market, x_min, x_max, y_min, y_max, width, height, margin)
    svg_elements.append(f"<circle cx='{x_market:.2f}' cy='{y_market:.2f}' r='6' fill='#08519c' stroke='black' stroke-width='1' />")
    svg_elements.append(
        f"<text x='{x_market + 10:.2f}' y='{y_market - 8:.2f}' font-size='14' fill='#08519c'>Market</text>"
    )

    svg_elements.append("</svg>")
    file_path.write_text("\n".join(svg_elements))


def main() -> None:
    data = read_market_data(DATA_FILE)
    risk_tolerance = 0.08
    risk_free_rate = 0.00015

    print("Question 1: mean vector (6 d.p.) and covariance matrix (6 d.p.)")
    for asset, mean in zip(data.assets, data.mean_vector):
        print(f"  {asset}: {mean:.6f}")
    print("Covariance matrix:")
    for row in data.covariance:
        print("  " + " ".join(f"{value:.6f}" for value in row))

    q2_weights, q2_mean, q2_std = constrained_optimal_weights(data.covariance, data.mean_vector, risk_tolerance)
    ones = [1.0] * len(data.assets)
    c_inv_r = gaussian_elimination_solve(data.covariance, data.mean_vector)
    c_inv_ones = gaussian_elimination_solve(data.covariance, ones)
    alpha = dot(ones, c_inv_r)
    beta = dot(ones, c_inv_ones)
    slope_terms = [c_inv_r[i] - (alpha / beta) * c_inv_ones[i] for i in range(len(data.assets))]
    intercept_terms = [c_inv_ones[i] / beta for i in range(len(data.assets))]
    q2_allocation = {asset: weight for asset, weight in zip(data.assets, q2_weights)}
    print("\nQuestion 2: optimal weights as functions of t (reported at t=0.08)")
    for asset, weight in q2_allocation.items():
        print(f"  {asset}: {weight:.6f}")
    print(f"  Expected return: {q2_mean:.6f}, risk (std): {q2_std:.6f}")
    print("  Weight(t) = slope * t + intercept")
    for asset, slope, intercept in zip(data.assets, slope_terms, intercept_terms):
        print(f"    {asset}: slope={slope:.6f}, intercept={intercept:.6f}")

    total_capital = 1_000_000.0
    q3_amounts = {asset: weight * total_capital for asset, weight in zip(data.assets, q2_weights)}
    print("\nQuestion 3: allocation for $1,000,000 (t=0.08)")
    for asset, amount in q3_amounts.items():
        print(f"  {asset}: ${amount:,.2f}")
    print(f"  Portfolio mean return: ${q2_mean * total_capital:,.2f}")
    print(f"  Portfolio risk (std): ${q2_std * total_capital:,.2f}")

    gmv_weights = global_minimum_variance_weights(data.covariance)
    gmv_mean, gmv_var, gmv_std = portfolio_statistics(gmv_weights, data.mean_vector, data.covariance)
    frontier = efficient_frontier(data.covariance, data.mean_vector)
    ensure_directory(FIGURE_DIR)
    generate_question4_svg(data, risk_tolerance, q2_weights, gmv_weights, frontier, FIGURE_DIR / "question4.svg")
    print("\nQuestion 4: figure saved to figures/question4.svg")
    print(f"  GMV mean: {gmv_mean:.6f}, GMV risk: {gmv_std:.6f}")

    q5_weights, q5_mean, q5_std = enumerate_no_short_solutions(data.covariance, data.mean_vector, risk_tolerance)
    print("\nQuestion 5a: optimal weights without short selling")
    for asset, weight in zip(data.assets, q5_weights):
        print(f"  {asset}: {weight:.6f}")
    print(f"  Expected return: {q5_mean:.6f}, risk (std): {q5_std:.6f}")

    lam, gamma = kkt_residuals(data.covariance, data.mean_vector, q5_weights, risk_tolerance)
    print("Question 5b: KKT verification")
    print(f"  Lambda: {lam:.6e}")
    for asset, value in zip(data.assets, gamma):
        print(f"  gamma_{asset}: {value:.6e}")

    risk_change = "lower" if q5_std < q2_std else "higher"
    return_change = "lower" if q5_mean < q2_mean else "higher"
    print(
        "Question 5c: Compared with Question 3, the constrained portfolio has "
        f"{risk_change} risk and {return_change} expected return."
    )

    tangent_weights, market_mean, market_std = tangent_portfolio(data.covariance, data.mean_vector, risk_free_rate)
    q7_risky_weights, q7_risk_free_weight, q7_mean, q7_std = risk_free_optimal_portfolio(
        data.covariance, data.mean_vector, risk_tolerance, risk_free_rate
    )
    generate_question7_svg(
        data,
        frontier,
        q2_weights,
        q7_risky_weights,
        q7_risk_free_weight,
        tangent_weights,
        risk_free_rate,
        (market_mean, q7_mean, market_std),
        FIGURE_DIR / "question7.svg",
    )

    print("\nQuestion 7a: optimal weights with risk-free asset")
    for asset, weight in zip(data.assets, q7_risky_weights):
        print(f"  {asset}: {weight:.6f}")
    print(f"  Risk-free asset weight: {q7_risk_free_weight:.6f}")
    print(f"  Expected return: {q7_mean:.6f}, risk (std): {q7_std:.6f}")

    print("Question 7b: figure saved to figures/question7.svg")
    print("  Market (tangency) portfolio and capital market line are included in the plot.")


if __name__ == "__main__":
    main()
