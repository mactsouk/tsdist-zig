# Time Series Distances Zig Library

A Zig library for calculating distances between time series, including Euclidean, Manhattan, Chebyshev, LCSS, EDR and Dynamic Time Warping (DTW) distances.

## Dynamic Time Warping (DTW) Performance

The DTW implementation includes an optional window parameter that constrains how far the warping path can deviate from the diagonal. This has significant performance implications:

### Without Window Constraint (`window = null`)
- Time Complexity: O(n*m)
- Space Complexity: O(n*m)
- Behavior: Examines all possible alignments between the two series
- Use Case: Only recommended for short time series (<100 points)

### With Window Constraint (`window = k`)
- Time Complexity: O(n*k)
- Space Complexity: O(n*k)
- Behavior: Only considers alignments within k steps of the diagonal
- Use Case: Recommended for longer series, provides a good balance between accuracy and performance

### Default Window
The `dtwDistanceDefaultWindow` function uses a window size of 10% of the maximum series length, which works well for most applications.

## Choosing a Window Size

1. **Small window (1-5)**: Use for nearly aligned series with little time shift
2. **Medium window (5-10% of series length)**: Good default for most applications
3. **Large window (>10% of series length)**: Use when significant time warping is expected
4. **No window (null)**: Only for short series or when exact DTW is required

## Memory Usage Warning

For long time series, the full DTW (without window) can consume significant memory:
- Two 1000-point series: ~8MB memory usage
- Two 10,000-point series: ~800MB memory usage
- Two 100,000-point series: ~80GB memory usage (likely to fail on most systems)

Always use a window constraint for series longer than a few hundred points.
