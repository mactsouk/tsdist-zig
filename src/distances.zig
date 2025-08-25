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
