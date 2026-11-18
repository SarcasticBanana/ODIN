function report=validateFinalPowerPlots(repositoryRoot,outputFolder)
% VALIDATEFINALPOWERPLOTS Native figure-preservation and graphics smoke test.
% report=validateFinalPowerPlots(repositoryRoot,outputFolder) uses the real
% ODIN cfg/workbook and one orbit, restores path/pwd, and leaves ALL figures
% open. outputFolder is an explicit artifact destination, never a default
% repository output directory. No protected source/workbook is modified.
if nargin<2
    error('finalPowerEPS:plotValidation:inputs','Supply repositoryRoot and an explicit outputFolder.');
end
repositoryRoot=char(repositoryRoot); outputFolder=char(outputFolder);
validateattributes(repositoryRoot,{'char'},{'row','nonempty'});
validateattributes(outputFolder,{'char'},{'row','nonempty'});
if ~isfolder(repositoryRoot) || ~isfile(fullfile(repositoryRoot,'odinPath.m'))
    error('finalPowerEPS:plotValidation:repository','repositoryRoot must be the ODIN-Model folder.');
end
[ok,repositoryAttribute]=fileattrib(repositoryRoot);
if ~ok, error('finalPowerEPS:plotValidation:repository','Could not resolve repositoryRoot.'); end
repositoryRoot=repositoryAttribute.Name;
originalPath=path; originalFolder=pwd;
cleanup=onCleanup(@()restoreSession(originalPath,originalFolder)); %#ok<NASGU>
packageRoot=fileparts(fileparts(mfilename('fullpath')));
% Resolve a relative artifact directory before changing the working folder.
if ~(startsWith(outputFolder,filesep) || ~isempty(regexp(outputFolder,'^[A-Za-z]:[\\/]','once')))
    outputFolder=fullfile(originalFolder,outputFolder);
end
cd(repositoryRoot);
odinPath;
addpath(packageRoot);
resolvedConfig=which('config'); resolvedPropagate=which('propagate');
expectedConfig=fullfile(repositoryRoot,'model','config.m');
expectedPropagate=fullfile(repositoryRoot,'model','propagate.m');
assert(strcmpi(resolvedConfig,expectedConfig),'Configuration did not resolve to the protected ODIN function.');
assert(strcmpi(resolvedPropagate,expectedPropagate),'Propagation did not resolve to the protected ODIN function.');
cfg=config(); out=propagate(cfg,1);
cfgBefore=cfg; outBefore=out;
existingFigures=findall(groot,'Type','figure');
options=struct('makePlots',true,'showSummary',false,'includePowerBalance',false, ...
    'arrayOverrides',struct('runIncidenceSweep',true, ...
    'incidenceAngles_deg',[0 45 90],'incidenceSummaryAngles_deg',[0 45 90]));
result=runFinalPowerEPS(cfg,out,options);
assert(isequaln(cfg,cfgBefore) && isequaln(out,outBefore),'Plot generation changed input data.');
assert(all(isgraphics(existingFigures,'figure')),'Package closed an existing figure.');
allFigures=findall(groot,'Type','figure');
newFigures=setdiff(allFigures,existingFigures);
assert(numel(newFigures)>=7,'Expected array, incidence, I-V and two diagnostic figures.');
assert(numel(result.figures)==2 && all(isgraphics(result.figures,'figure')), ...
    'Both added diagnostic figures must be available.');
assert(all(ismember(result.figures,newFigures)),'Diagnostic handles must refer to newly created figures.');
assert(strcmpi(which('config'),resolvedConfig) && strcmpi(which('propagate'),resolvedPropagate), ...
    'Package changed shared-function resolution.');
if ~isfolder(outputFolder), mkdir(outputFolder); end
temperaturePNG=fullfile(outputFolder,'final-power-eps-temperature-difference.png');
operationPNG=fullfile(outputFolder,'final-power-eps-operation.png');
exportgraphics(result.figures(1),temperaturePNG,'Resolution',120);
exportgraphics(result.figures(2),operationPNG,'Resolution',120);
assert(isfile(temperaturePNG) && isfile(operationPNG),'Native PNG export did not complete.');
names=cell(numel(newFigures),1);
for k=1:numel(newFigures), names{k}=newFigures(k).Name; end
report=struct('passed',true,'originalFigureCount',numel(existingFigures), ...
    'newFigureCount',numel(newFigures),'newFigureNames',{names}, ...
    'temperaturePNG',temperaturePNG,'operationPNG',operationPNG, ...
    'sharedConfig',resolvedConfig,'sharedPropagate',resolvedPropagate);
fprintf('PASS: %d existing figures retained; %d package figures rendered; two native diagnostic exports created.\n', ...
    report.originalFigureCount,report.newFigureCount);
end
function restoreSession(originalPath,originalFolder)
path(originalPath);
cd(originalFolder);
end
