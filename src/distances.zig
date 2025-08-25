const std = @import("std");
const math = std.math;

pub const Error = error{InvalidInput};

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

// Space-optimized DTW using only two rows - always returns finite results
pub fn dtwDistance(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
    window: ?usize,
) (std.mem.Allocator.Error || Error)!f64 {
    const n = a.len;
    const m = b.len;

    if (n == 0 or m == 0) return Error.InvalidInput;

    const cols = m + 1;

    // Only need two rows: previous and current
    var prev_row = try allocator.alloc(f64, cols);
    defer allocator.free(prev_row);
    var curr_row = try allocator.alloc(f64, cols);
    defer allocator.free(curr_row);

    // Initialize first row with cumulative costs (aligning empty sequence with b)
    prev_row[0] = 0.0;
    for (1..cols) |j| {
        prev_row[j] = prev_row[j - 1] + @abs(b[j - 1]);
    }

    // Process each row
    for (1..(n + 1)) |i| {
        // First column: cumulative cost of aligning a with empty sequence
        curr_row[0] = prev_row[0] + @abs(a[i - 1]);

        for (1..cols) |j| {
            // Apply window constraint if specified
            if (window) |w| {
                const i_int: i64 = @intCast(i);
                const j_int: i64 = @intCast(j);
                if (@abs(i_int - j_int) > @as(i64, @intCast(w))) {
                    // Outside window - use point-to-point distance as fallback
                    curr_row[j] = @abs(a[i - 1] - b[j - 1]);
                    continue;
                }
            }

            // DTW recurrence relation
            const cost = @abs(a[i - 1] - b[j - 1]);
            const deletion = prev_row[j]; // from above
            const insertion = curr_row[j - 1]; // from left
            const match = prev_row[j - 1]; // from diagonal

            const min_val = @min(@min(deletion, insertion), match);
            curr_row[j] = cost + min_val;
        }

        // Swap rows for next iteration
        const temp = prev_row;
        prev_row = curr_row;
        curr_row = temp;
    }

    return prev_row[m];
}

// Convenience wrappers
pub fn dtwDistanceDefaultWindow(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
) (std.mem.Allocator.Error || Error)!f64 {
    const maxlen = if (a.len > b.len) a.len else b.len;
    const minlen = if (a.len < b.len) a.len else b.len;

    // Ensure window is at least as large as the length difference
    const length_diff = maxlen - minlen;
    const default_window = maxlen / 10;
    const window_size = @max(@max(default_window, length_diff), 1);

    return dtwDistance(allocator, a, b, window_size);
}

pub fn dtwDistanceFull(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
) (std.mem.Allocator.Error || Error)!f64 {
    return dtwDistance(allocator, a, b, null);
}

// Space-optimized LCSS using only two rows
pub fn lcssDistance(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
    epsilon: f64,
) (std.mem.Allocator.Error || Error)!f64 {
    const n = a.len;
    const m = b.len;

    if (n == 0 and m == 0) return 0.0;
    if (n == 0 or m == 0) return 1.0;

    const cols = m + 1;

    var prev_row = try allocator.alloc(u32, cols);
    defer allocator.free(prev_row);
    var curr_row = try allocator.alloc(u32, cols);
    defer allocator.free(curr_row);

    // Initialize first row to 0
    for (0..cols) |j| prev_row[j] = 0;

    // Process each row
    for (1..(n + 1)) |i| {
        curr_row[0] = 0; // First column is always 0

        for (1..cols) |j| {
            if (@abs(a[i - 1] - b[j - 1]) <= epsilon) {
                curr_row[j] = prev_row[j - 1] + 1;
            } else {
                const left = curr_row[j - 1];
                const top = prev_row[j];
                curr_row[j] = if (left > top) left else top;
            }
        }

        // Swap rows for next iteration
        const temp = prev_row;
        prev_row = curr_row;
        curr_row = temp;
    }

    const lcs_len: f64 = @floatFromInt(prev_row[m]);
    const max_len: f64 = @floatFromInt(if (n > m) n else m);
    return 1.0 - (lcs_len / max_len);
}

// Fast LCSS for cases where you only need approximate results
pub fn lcssDistanceFast(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
    epsilon: f64,
    max_length: ?usize, // Limit processing to first max_length elements
) (std.mem.Allocator.Error || Error)!f64 {
    const limit = if (max_length) |lim| lim else @min(a.len, b.len);
    const n = @min(a.len, limit);
    const m = @min(b.len, limit);

    const a_slice = a[0..n];
    const b_slice = b[0..m];

    return lcssDistance(allocator, a_slice, b_slice, epsilon);
}

// Space-optimized EDR using only two rows
pub fn edrDistance(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
    epsilon: f64,
) (std.mem.Allocator.Error || Error)!f64 {
    const n = a.len;
    const m = b.len;

    if (n == 0) return @floatFromInt(m);
    if (m == 0) return @floatFromInt(n);

    const cols = m + 1;

    var prev_row = try allocator.alloc(f32, cols); // Use f32 for memory efficiency
    defer allocator.free(prev_row);
    var curr_row = try allocator.alloc(f32, cols);
    defer allocator.free(curr_row);

    // Initialize first row: [0, 1, 2, 3, ...]
    for (0..cols) |j| prev_row[j] = @floatFromInt(j);

    // Process each row
    for (1..(n + 1)) |i| {
        curr_row[0] = @floatFromInt(i); // First column: [0, 1, 2, 3, ...]

        for (1..cols) |j| {
            const cost: f32 = if (@abs(a[i - 1] - b[j - 1]) <= epsilon) 0.0 else 1.0;
            const del = prev_row[j] + 1.0; // deletion
            const ins = curr_row[j - 1] + 1.0; // insertion
            const sub = prev_row[j - 1] + cost; // substitution
            curr_row[j] = @min(@min(del, ins), sub);
        }

        // Swap rows for next iteration
        const temp = prev_row;
        prev_row = curr_row;
        curr_row = temp;
    }

    return @floatCast(prev_row[m]);
}

// Fast EDR for large sequences - limits processing to first N elements
pub fn edrDistanceFast(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
    epsilon: f64,
    max_length: ?usize,
) (std.mem.Allocator.Error || Error)!f64 {
    const limit = if (max_length) |lim| lim else @min(a.len, b.len);
    const n = @min(a.len, limit);
    const m = @min(b.len, limit);

    const a_slice = a[0..n];
    const b_slice = b[0..m];

    return edrDistance(allocator, a_slice, b_slice, epsilon);
}

// Very fast EDR using sampling for huge sequences
pub fn edrDistanceSampled(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
    epsilon: f64,
    sample_rate: usize, // Process every Nth element
) (std.mem.Allocator.Error || Error)!f64 {
    if (sample_rate == 0) return Error.InvalidInput;

    const n_sampled = (a.len + sample_rate - 1) / sample_rate; // Ceiling division
    const m_sampled = (b.len + sample_rate - 1) / sample_rate;

    var a_sampled = try allocator.alloc(f64, n_sampled);
    defer allocator.free(a_sampled);
    var b_sampled = try allocator.alloc(f64, m_sampled);
    defer allocator.free(b_sampled);

    // Sample elements
    for (0..n_sampled) |i| {
        const sample_idx = i * sample_rate;
        a_sampled[i] = if (sample_idx < a.len) a[sample_idx] else a[a.len - 1];
    }

    for (0..m_sampled) |i| {
        const sample_idx = i * sample_rate;
        b_sampled[i] = if (sample_idx < b.len) b[sample_idx] else b[b.len - 1];
    }

    return edrDistance(allocator, a_sampled, b_sampled, epsilon);
}
