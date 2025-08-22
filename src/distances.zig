const std = @import("std");
const math = std.math;

pub const Error = error{InvalidInput};

// Helper: compute index into flat (rows x cols) matrix
fn idx(cols: usize, i: usize, j: usize) usize {
    return i * cols + j;
}

// Equal-length distance functions
pub fn euclideanDistance(a: []const f64, b: []const f64) !f64 {
    if (a.len != b.len) return Error.InvalidInput;
    var sum: f64 = 0.0;
    for (0..a.len) |i| {
        const diff = a[i] - b[i];
        sum += diff * diff;
    }
    return math.sqrt(sum);
}

pub fn manhattanDistance(a: []const f64, b: []const f64) !f64 {
    if (a.len != b.len) return Error.InvalidInput;
    var sum: f64 = 0.0;
    for (0..a.len) |i| {
        sum += @abs(a[i] - b[i]);
    }
    return sum;
}

pub fn chebyshevDistance(a: []const f64, b: []const f64) !f64 {
    if (a.len != b.len) return Error.InvalidInput;
    var max_diff: f64 = 0.0;
    for (0..a.len) |i| {
        const diff = @abs(a[i] - b[i]);
        if (diff > max_diff) max_diff = diff;
    }
    return max_diff;
}

// Dynamic Time Warping with window constraint
pub fn dtwDistance(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
    window: ?usize,
) !f64 {
    const n = a.len;
    const m = b.len;

    // Determine the window size
    const w: usize = if (window) |win| win else if (n > m) n else m;

    const rows = n + 1;
    const cols = m + 1;
    const total = rows * cols;

    var mat = try allocator.alloc(f64, total);
    defer allocator.free(mat);

    // initialize to +inf
    for (0..total) |i| {
        mat[i] = math.inf(f64);
    }

    // base
    mat[idx(cols, 0, 0)] = 0.0;

    for (1..rows) |i| {
        const j_start = if (i > w) i - w else 1;
        const j_end = if ((i + w) < m) i + w else m;
        if (j_start > m) continue;

        for (j_start..(j_end + 1)) |j| {
            const cost = @abs(a[i - 1] - b[j - 1]);
            const v1 = mat[idx(cols, i - 1, j)];
            const v2 = mat[idx(cols, i, j - 1)];
            const v3 = mat[idx(cols, i - 1, j - 1)];
            const min_val = @min(@min(v1, v2), v3);
            mat[idx(cols, i, j)] = cost + min_val;
        }
    }

    return mat[idx(cols, n, m)];
}

// Convenience wrappers
pub fn dtwDistanceDefaultWindow(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
) !f64 {
    const maxlen = if (a.len > b.len) a.len else b.len;
    const window_size = if (maxlen / 10 > 1) maxlen / 10 else 1;
    return dtwDistance(allocator, a, b, window_size);
}

pub fn dtwDistanceFull(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
) !f64 {
    return dtwDistance(allocator, a, b, null);
}

// Longest Common Subsequence (LCSS) distance
pub fn lcssDistance(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
    epsilon: f64,
) !f64 {
    const n = a.len;
    const m = b.len;
    const rows = n + 1;
    const cols = m + 1;
    const total = rows * cols;

    var mat = try allocator.alloc(usize, total);
    defer allocator.free(mat);

    // initialize to 0
    for (0..total) |i| mat[i] = 0;

    for (1..rows) |i| {
        for (1..cols) |j| {
            if (@abs(a[i - 1] - b[j - 1]) <= epsilon) {
                mat[idx(cols, i, j)] = mat[idx(cols, i - 1, j - 1)] + 1;
            } else {
                const left = mat[idx(cols, i, j - 1)];
                const top = mat[idx(cols, i - 1, j)];
                mat[idx(cols, i, j)] = if (left > top) left else top;
            }
        }
    }

    const lcs_len: f64 = @floatFromInt(mat[idx(cols, n, m)]);
    const max_len: f64 = @floatFromInt(if (n > m) n else m);
    if (max_len == 0.0) return 0.0;
    return 1.0 - (lcs_len / max_len);
}

// Edit Distance on Real sequences (EDR)
pub fn edrDistance(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
    epsilon: f64,
) !f64 {
    const n = a.len;
    const m = b.len;
    const rows = n + 1;
    const cols = m + 1;
    const total = rows * cols;

    var mat = try allocator.alloc(f64, total);
    defer allocator.free(mat);

    // initialize
    for (0..total) |i| mat[i] = 0.0;

    // first row and first column
    for (0..cols) |j| mat[idx(cols, 0, j)] = @floatFromInt(j);
    for (1..rows) |i| mat[idx(cols, i, 0)] = @floatFromInt(i);

    for (1..rows) |i| {
        for (1..cols) |j| {
            const cost: f64 = if (@abs(a[i - 1] - b[j - 1]) <= epsilon) 0 else 1;
            const del = mat[idx(cols, i - 1, j)] + 1.0;
            const ins = mat[idx(cols, i, j - 1)] + 1.0;
            const sub = mat[idx(cols, i - 1, j - 1)] + cost;
            mat[idx(cols, i, j)] = @min(@min(del, ins), sub);
        }
    }

    return mat[idx(cols, n, m)];
}
