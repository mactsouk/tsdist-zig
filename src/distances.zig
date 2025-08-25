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

// DTW initialization strategies
pub const DTWInit = enum {
    Infinity,           // Traditional: initialize with infinity
    Cumulative,         // Initialize with cumulative distances
    Euclidean,          // Initialize with point-to-point Euclidean distances
    Manhattan,          // Initialize with point-to-point Manhattan distances
};

// Fixed Dynamic Time Warping with window constraint and initialization options
pub fn dtwDistanceWithInit(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
    window: ?usize,
    init_strategy: DTWInit,
) (std.mem.Allocator.Error || Error)!f64 {
    const n = a.len;
    const m = b.len;
    
    if (n == 0 or m == 0) return Error.InvalidInput;

    const rows = n + 1;
    const cols = m + 1;
    const total = rows * cols;

    var mat = try allocator.alloc(f64, total);
    defer allocator.free(mat);

    // Initialize matrix based on strategy
    switch (init_strategy) {
        .Infinity => {
            // Traditional initialization
            for (0..total) |i| {
                mat[i] = math.inf(f64);
            }
            mat[idx(cols, 0, 0)] = 0.0;
        },
        .Cumulative => {
            // Initialize with cumulative distances
            mat[idx(cols, 0, 0)] = 0.0;
            
            // First row: cumulative cost of aligning empty sequence with b
            for (1..cols) |j| {
                mat[idx(cols, 0, j)] = mat[idx(cols, 0, j - 1)] + @abs(b[j - 1]);
            }
            
            // First column: cumulative cost of aligning a with empty sequence  
            for (1..rows) |i| {
                mat[idx(cols, i, 0)] = mat[idx(cols, i - 1, 0)] + @abs(a[i - 1]);
            }
            
            // Initialize remaining cells to infinity
            for (1..rows) |i| {
                for (1..cols) |j| {
                    mat[idx(cols, i, j)] = math.inf(f64);
                }
            }
        },
        .Euclidean => {
            // Initialize with point-to-point Euclidean distances
            mat[idx(cols, 0, 0)] = 0.0;
            
            // First row and column as cumulative
            for (1..cols) |j| {
                mat[idx(cols, 0, j)] = mat[idx(cols, 0, j - 1)] + @abs(b[j - 1]);
            }
            for (1..rows) |i| {
                mat[idx(cols, i, 0)] = mat[idx(cols, i - 1, 0)] + @abs(a[i - 1]);
            }
            
            // Initialize interior with Euclidean distances from origin
            for (1..rows) |i| {
                for (1..cols) |j| {
                    const dist_a = @abs(a[i - 1]);
                    const dist_b = @abs(b[j - 1]);
                    mat[idx(cols, i, j)] = math.sqrt(dist_a * dist_a + dist_b * dist_b);
                }
            }
        },
        .Manhattan => {
            // Initialize with Manhattan distances
            mat[idx(cols, 0, 0)] = 0.0;
            
            for (1..cols) |j| {
                mat[idx(cols, 0, j)] = mat[idx(cols, 0, j - 1)] + @abs(b[j - 1]);
            }
            for (1..rows) |i| {
                mat[idx(cols, i, 0)] = mat[idx(cols, i - 1, 0)] + @abs(a[i - 1]);
            }
            
            for (1..rows) |i| {
                for (1..cols) |j| {
                    mat[idx(cols, i, j)] = @abs(a[i - 1]) + @abs(b[j - 1]);
                }
            }
        },
    }

    // Rest of DTW algorithm remains the same...

    // If no window constraint, use full DTW
    if (window == null) {
        for (1..rows) |i| {
            for (1..cols) |j| {
                const cost = @abs(a[i - 1] - b[j - 1]);
                const v1 = mat[idx(cols, i - 1, j)];
                const v2 = mat[idx(cols, i, j - 1)];
                const v3 = mat[idx(cols, i - 1, j - 1)];
                const min_val = @min(@min(v1, v2), v3);
                if (min_val != math.inf(f64)) {
                    mat[idx(cols, i, j)] = cost + min_val;
                }
            }
        }
        return mat[idx(cols, n, m)];
    }

    // Windowed DTW with proper Sakoe-Chiba band
    const w = window.?;
    
    for (1..rows) |i| {
        for (1..cols) |j| {
            // Sakoe-Chiba band constraint: |i - j| <= w
            const i_int: i64 = @intCast(i);
            const j_int: i64 = @intCast(j);
            const diff = @abs(i_int - j_int);
            const w_int: i64 = @intCast(w);
            
            if (diff <= w_int) {
                const cost = @abs(a[i - 1] - b[j - 1]);
                const v1 = mat[idx(cols, i - 1, j)];
                const v2 = mat[idx(cols, i, j - 1)];
                const v3 = mat[idx(cols, i - 1, j - 1)];
                const min_val = @min(@min(v1, v2), v3);
                
                if (min_val != math.inf(f64)) {
                    mat[idx(cols, i, j)] = cost + min_val;
                }
            }
        }
    }

    // Check if final cell is reachable
    const result = mat[idx(cols, n, m)];
    if (result == math.inf(f64)) {
        // Window too restrictive - fall back to larger window
        const fallback_window = @max(w * 2, @abs(@as(i64, @intCast(n)) - @as(i64, @intCast(m))));
        return dtwDistance(allocator, a, b, fallback_window);
    }
    
    return result;
}

// Original DTW function with default Euclidean initialization
pub fn dtwDistance(
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
    window: ?usize,
) (std.mem.Allocator.Error || Error)!f64 {
    return dtwDistanceWithInit(allocator, a, b, window, .Euclidean);
}

// Convenience wrappers with better defaults
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

// Longest Common Subsequence (LCSS) distance
// Optimized Longest Common Subsequence (LCSS) distance with space optimization
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

    // Space optimization: use only two rows instead of full matrix
    // This reduces space complexity from O(n*m) to O(min(n,m))
    const use_cols = if (m <= n) m else n;
    const use_rows = if (m <= n) n else m;
    const should_transpose = m > n;
    
    const cols = use_cols + 1;
    
    var prev_row = try allocator.alloc(u32, cols);
    defer allocator.free(prev_row);
    var curr_row = try allocator.alloc(u32, cols);
    defer allocator.free(curr_row);
    
    // Initialize first row to 0
    for (0..cols) |j| prev_row[j] = 0;
    
    // Process each row
    for (1..(use_rows + 1)) |i| {
        curr_row[0] = 0; // First column is always 0
        
        for (1..cols) |j| {
            // Handle transposition if needed
            const a_val = if (should_transpose) b[i - 1] else a[i - 1];
            const b_val = if (should_transpose) a[j - 1] else b[j - 1];
            
            if (@abs(a_val - b_val) <= epsilon) {
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
    
    const lcs_len: f64 = @floatFromInt(prev_row[use_cols]);
    const max_len: f64 = @floatFromInt(use_rows);
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

// Edit Distance on Real sequences (EDR)
// Optimized Edit Distance on Real sequences (EDR) with space optimization
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

    // Space optimization: use only two rows instead of full matrix
    // This reduces space complexity from O(n*m) to O(min(n,m))
    const use_cols = if (m <= n) m else n;
    const use_rows = if (m <= n) n else m;
    const should_transpose = m > n;
    
    const cols = use_cols + 1;
    
    var prev_row = try allocator.alloc(f32, cols); // Use f32 instead of f64 for memory efficiency
    defer allocator.free(prev_row);
    var curr_row = try allocator.alloc(f32, cols);
    defer allocator.free(curr_row);
    
    // Initialize first row: [0, 1, 2, 3, ...]
    for (0..cols) |j| prev_row[j] = @floatFromInt(j);
    
    // Process each row
    for (1..(use_rows + 1)) |i| {
        curr_row[0] = @floatFromInt(i); // First column: [0, 1, 2, 3, ...]
        
        for (1..cols) |j| {
            // Handle transposition if needed
            const a_val = if (should_transpose) b[i - 1] else a[i - 1];
            const b_val = if (should_transpose) a[j - 1] else b[j - 1];
            
            const cost: f32 = if (@abs(a_val - b_val) <= epsilon) 0.0 else 1.0;
            const del = prev_row[j] + 1.0;      // deletion
            const ins = curr_row[j - 1] + 1.0;  // insertion
            const sub = prev_row[j - 1] + cost; // substitution
            curr_row[j] = @min(@min(del, ins), sub);
        }
        
        // Swap rows for next iteration
        const temp = prev_row;
        prev_row = curr_row;
        curr_row = temp;
    }
    
    return @floatCast(prev_row[use_cols]);
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
