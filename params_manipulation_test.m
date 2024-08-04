function tests = params_manipulation_test
tests = functiontests(localfunctions);
end

function test_params2inds_1(testCase)
params = testCase.TestData.params;
prms.f1 = 1;
prms.f2.v1.f1 = ["a"];
prms.f2.v1.f2 = 9;
real_inds = [1 1 1 1];
computed_inds = params2inds(prms,params);
verifyEqual(testCase,computed_inds,real_inds)
end

function test_params2inds_2(testCase)
params = testCase.TestData.params;
prms.f1 = 5;
prms.f2.v2.f1 = ["cc"];
real_inds = [5 2 3];
computed_inds = params2inds(prms,params);
verifyEqual(testCase,computed_inds,real_inds)
end


function test_params2inds_invalidValue(testCase)
params = testCase.TestData.params;
prms.f1 = 6;
prms.f2.v1.f1 = ["a"];
prms.f2.v1.f2 = 9;
verifyError(testCase,@()params2inds(prms,params),'params2inds:invalid_value')
end

function test_params2inds_missingField(testCase)
params = testCase.TestData.params;
prms.f1 = 5;
prms.f2.v1.f1 = ["a"];
% prms.f2.v1.f2 = 9;
verifyError(testCase,@()params2inds(prms,params),'params2inds:missing_field')
end

function test_inds2params_1(testCase)
params = testCase.TestData.params;
real_prms.f1 = 1;
real_prms.f2.v1.f1 = ["a"];
real_prms.f2.v1.f2 = 9;
inds = [1 1 1 1];
computed_prms = inds2params(inds,params);
verifyEqual(testCase,computed_prms,real_prms)
end

function test_inds2params_2(testCase)
params = testCase.TestData.params;
real_prms.f1 = 5;
real_prms.f2.v2.f1 = ["cc"];
inds = [5 2 3];
computed_prms = inds2params(inds,params);
verifyEqual(testCase,computed_prms,real_prms)
end


function test_inds2params_outOfRange(testCase)
params = testCase.TestData.params;
inds = [5 3];
verifyError(testCase,@()inds2params(inds,params),'inds2params:outOfRange')
end

function test_inds2params_insufficientInds(testCase)
params = testCase.TestData.params;
inds = [5 2];
verifyError(testCase,@()inds2params(inds,params),'inds2params:insufficientInds')
end

function test_inds2params_tooManyInds(testCase)
params = testCase.TestData.params;
inds = [1 1 1 1 1];
verifyError(testCase,@()inds2params(inds,params),'inds2params:tooManyInds')
end

function test_generate_combinations(testCase)
params = testCase.TestData.params;
[params_all,inds_all] = generate_combinations(params);
verifyEqual(testCase,numel(params_all),numel(inds_all))
for i = 1:numel(params_all)
    verifyEqual(testCase,params2inds(params_all{i}.params,params),inds_all{i})
end
end

function setupOnce(testCase)  % do not change function name
params.f1 = 1:5;
params.f2.v1.f1 = ["a", "b"];
params.f2.v1.f2 = 9:20;
params.f2.v2.f1 = ["aa", "bb", "cc"];

testCase.TestData.params = params;
end