const std = @import("std");
const distances = @import("distances");

test "euclidean distance" {
    const a = [_]f64{ 1.0, 2.0, 3.0 };
    const b = [_]f64{ 4.0, 5.0, 6.0 };
    const result = try distances.euclideanDistance(&a, &b);
    try std.testing.expectApproxEqAbs(@sqrt(27.0), result, 1e-10);
}

test "manhattan distance" {
    const a = [_]f64{ 1.0, 2.0, 3.0 };
    const b = [_]f64{ 4.0, 5.0, 6.0 };
    const result = try distances.manhattanDistance(&a, &b);
    try std.testing.expectApproxEqAbs(9.0, result, 1e-10);
}

test "chebyshev distance" {
    const a = [_]f64{ 1.0, 2.0, 8.0 };
    const b = [_]f64{ 4.0, 5.0, 6.0 };
    const result = try distances.chebyshevDistance(&a, &b);
    try std.testing.expectApproxEqAbs(2.0, result, 1e-10);
}

test "invalid input" {
    const a = [_]f64{ 1.0, 2.0 };
    const b = [_]f64{1.0};
    try std.testing.expectError(error.InvalidInput, distances.euclideanDistance(&a, &b));
    try std.testing.expectError(error.InvalidInput, distances.manhattanDistance(&a, &b));
    try std.testing.expectError(error.InvalidInput, distances.chebyshevDistance(&a, &b));
}

test "dtw distance" {
    const allocator = std.testing.allocator;
    const a = [_]f64{ 1.0, 2.0, 3.0 };
    const b = [_]f64{ 2.0, 3.0, 4.0, 5.0 };
    const result = try distances.dtwDistance(allocator, &a, &b, null);
    try std.testing.expectApproxEqAbs(@sqrt(7.0), result, 1e-10);
}

test "dtw distance with window" {
    const allocator = std.testing.allocator;
    const a = [_]f64{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    const b = [_]f64{ 2.0, 3.0, 4.0, 5.0, 6.0 };

    const result_full = try distances.dtwDistance(allocator, &a, &b, null);
    const result_window = try distances.dtwDistance(allocator, &a, &b, 2);

    try std.testing.expectApproxEqAbs(result_full, result_window, 1e-10);
}

test "lcss distance" {
    const allocator = std.testing.allocator;
    const a = [_]f64{ 1.0, 2.0, 3.0, 4.0 };
    const b = [_]f64{ 2.0, 3.0, 5.0 };
    const result = try distances.lcssDistance(allocator, &a, &b, 0.5);
    try std.testing.expectApproxEqAbs(0.5, result, 1e-10);
}

test "edr distance" {
    const allocator = std.testing.allocator;
    const a = [_]f64{ 1.0, 2.0, 3.0 };
    const b = [_]f64{ 2.0, 3.0, 4.0, 5.0 };
    const result = try distances.edrDistance(allocator, &a, &b, 0.5);
    try std.testing.expectApproxEqAbs(3.0, result, 1e-10);
}
