classdef BassTestCase < matlab.unittest.TestCase
    % Shared base class for BASS tests. Puts the project root (the parent of
    % this tests/ folder, which contains the source .m files) on the MATLAB
    % path for the duration of each test class, so tests pass regardless of
    % the directory they are launched from.

    methods (TestClassSetup)
        function addProjectRootToPath(tc)
            import matlab.unittest.fixtures.PathFixture
            thisDir = fileparts(mfilename('fullpath'));
            projectRoot = fileparts(thisDir);
            tc.applyFixture(PathFixture(projectRoot));
        end
    end
end
