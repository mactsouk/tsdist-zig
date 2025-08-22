const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Library module
    const distances_module = b.addModule("distances", .{
        .root_source_file = b.path("src/distances.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Library artifact (static)
    const lib = b.addLibrary(.{
        .linkage = .static,
        .name = "distances",
        .root_module = distances_module,
    });
    b.installArtifact(lib);

    // Test module — points at tests/ (plural)
    const test_root = b.createModule(.{
        .root_source_file = b.path("tests/distances_test.zig"),
        .target = target,
        .optimize = optimize,
    });

    const tests = b.addTest(.{ .root_module = test_root });
    // Allow tests to import the library as `@import("distances")`
    tests.root_module.addImport("distances", distances_module);

    const run_tests = b.addRunArtifact(tests);
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_tests.step);
}
