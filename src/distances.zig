const std = @import("std");
const math = std.math;

// Equal-length distance functions (unchanged)
pub fn euclideanDistance(a: []const f64, b: []const f64) !f64 {
    if (a.len != b.len) return error.InvalidInput;
    var sum: f64 = 0.0;
    for (a, b) |x, y| {
        const diff = x - y;
        sum += diff * diff;
    }
    return math.sqrt(sum);
}

pub fn manhattanDistance(a: []const f64, b: []const f64) !f64 {
    if (a.len != b.len) return error.InvalidInput;
    var sum: f64 = 0.0;
    for (a, b) |x, y| {
        sum += @abs(x - y);
    }
    return sum;
}

pub fn chebyshevDistance(a: []const f64, b: []const f64) !f64 {
    if (a.len != b.len) return error.InvalidInput;
    var max_diff: f64 = 0.0;
    for (a, b) |x, y| {
        const diff = @abs(x - y);
        if (diff > max_diff) max_diff = diff;
    }
    return max_diff;
}

// Dynamic Time Warping with window constraint
pub fn dtwDistance(allocator: std.mem.Allocator, a: []const f64, b: []const f64, window: ?usize) !f64 {
    const n = a.len;
    const m = b.len;

    // Determine the window size
    const w = if (window) |win| win else math.max(n, m);

    // Create a matrix for dynamic programming
    var matrix = try allocator.alloc([]f64, n + 1);
    // Ensure proper cleanup even on errors
    errdefer {
        for (matrix[0..]) |row| {
            allocator.free(row);
        }
        allocator.free(matrix);
    }

    for (0..n + 1) |i| {
        matrix[i] = try allocator.alloc(f64, m + 1);
        @memset(matrix[i], std.math.inf(f64));

        // Ensure proper cleanup for rows on errors
        errdefer {
            for (matrix[0..i]) |row| {
                allocator.free(row);
            }
            allocator.free(matrix);
        }
    }

    // Final cleanup
    defer {
        for (matrix) |row| {
            allocator.free(row);
        }
        allocator.free(matrix);
    }

    matrix[0][0] = 0;

    for (1..n + 1) |i| {
        // Calculate the window boundaries for this row
        const j_start = if (i > w) i - w else 1;
        const j_end = math.min(m, i + w) + 1;

        // Skip if window is outside bounds
        if (j_start > m) continue;

        for (j_start..j_end) |j| {
            const cost = @abs(a[i - 1] - b[j - 1]);
            const min_val = math.min(math.min(matrix[i - 1][j], matrix[i][j - 1]), matrix[i - 1][j - 1]);
            matrix[i][j] = cost + min_val;
        }
    }

    return matrix[n][m];
}

// Convenience functions for DTW
pub fn dtwDistanceDefaultWindow(allocator: std.mem.Allocator, a: []const f64, b: []const f64) !f64 {
    const window_size = math.max(1, math.max(a.len, b.len) / 10);
    return dtwDistance(allocator, a, b, window_size);
}

pub fn dtwDistanceFull(allocator: std.mem.Allocator, a: []const f64, b: []const f64) !f64 {
    return dtwDistance(allocator, a, b, null);
}

// Longest Common Subsequence distance
pub fn lcssDistance(allocator: std.mem.Allocator, a: []const f64, b: []const f64, epsilon: f64) !f64 {
    const n = a.len;
    const m = b.len;

    // Create a matrix for dynamic programming
    var matrix = try allocator.alloc([]usize, n + 1);
    errdefer {
        for (matrix[0..]) |row| allocator.free(row);
        allocator.free(matrix);
    }

    for (0..n + 1) |i| {
        matrix[i] = try allocator.alloc(usize, m + 1);
        @memset(matrix[i], 0);

        errdefer {
            for (matrix[0..i]) |row| allocator.free(row);
            allocator.free(matrix);
        }
    }

    defer {
        for (matrix) |row| allocator.free(row);
        allocator.free(matrix);
    }

    for (1..n + 1) |i| {
        for (1..m + 1) |j| {
            if (@abs(a[i - 1] - b[j - 1]) <= epsilon) {
                matrix[i][j] = matrix[i - 1][j - 1] + 1;
            } else {
                matrix[i][j] = math.max(matrix[i][j - 1], matrix[i - 1][j]);
            }
        }
    }

    const lcs_length = @as(f64, @floatFromInt(matrix[n][m]));
    const max_len = @as(f64, @floatFromInt(math.max(n, m)));
    return 1.0 - (lcs_length / max_len);
}

// Edit Distance on Real sequence
pub fn edrDistance(allocator: std.mem.Allocator, a: []const f64, b: []const f64, epsilon: f64) !f64 {
    const n = a.len;
    const m = b.len;

    // Create a matrix for dynamic programming
    var matrix = try allocator.alloc([]f64, n + 1);
    errdefer {
        for (matrix[0..]) |row| allocator.free(row);
        allocator.free(matrix);
    }

    for (0..n + 1) |i| {
        matrix[i] = try allocator.alloc(f64, m + 1);

        if (i == 0) {
            for (0..m + 1) |j| {
                matrix[i][j] = @as(f64, @floatFromInt(j));
            }
        } else {
            matrix[i][0] = @as(f64, @floatFromInt(i));
            // Initialize the rest to 0
            for (1..m + 1) |j| {
                matrix[i][j] = 0;
            }
        }

        errdefer {
            for (matrix[0..i]) |row| allocator.free(row);
            allocator.free(matrix);
        }
    }

    defer {
        for (matrix) |row| allocator.free(row);
        allocator.free(matrix);
    }

    for (1..n + 1) |i| {
        for (1..m + 1) |j| {
            const cost = if (@abs(a[i - 1] - b[j - 1]) <= epsilon) 0 else 1;
            matrix[i][j] = math.min(math.min(matrix[i - 1][j] + 1, matrix[i][j - 1] + 1), matrix[i - 1][j - 1] + cost);
        }
    }

    return matrix[n][m];
}
