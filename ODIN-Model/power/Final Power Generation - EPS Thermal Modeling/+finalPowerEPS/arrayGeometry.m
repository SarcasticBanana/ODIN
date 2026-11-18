function g=arrayGeometry(arrayArea_m2,cellCoverage)
% ARRAYGEOMETRY Shared provisional geometry for the NanoAvionics 12U model.
% Areas refer to ONE collecting front per panel, summed over the array.
% Cell count is user supplied. New cell dimensions/active area are NOT vendor verified.
% Default preserves the coverage inferred from the OLD CAD areas under uniform resizing.
g=finalPowerEPS.modelSettings('geometry');
if nargin<1 || isempty(arrayArea_m2), arrayArea_m2=g.defaultArrayArea_m2; end
if nargin<2 || isempty(cellCoverage)
    cellCoverage=g.totalCells*g.legacyCellArea_m2/g.legacyArrayArea_m2;
end
validateattributes(arrayArea_m2,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(cellCoverage,{'numeric'},{'real','finite','scalar','>',0,'<=',1});
g.arrayArea_m2=arrayArea_m2; g.cellCoverage=cellCoverage;
g.totalCellArea_m2=arrayArea_m2*cellCoverage;
g.cellArea_m2=g.totalCellArea_m2/g.totalCells;
g.cellArea_cm2=g.cellArea_m2*1e4;
g.maximumCellArea_m2=arrayArea_m2/g.totalCells;
g.uncoveredArea_m2=arrayArea_m2-g.totalCellArea_m2;
g.oldCellsTotalArea_m2=g.totalCells*g.legacyCellArea_m2;
g.oldCellsFit=g.oldCellsTotalArea_m2<=arrayArea_m2;
g.uniformLinearScale=sqrt(g.cellArea_m2/g.legacyCellArea_m2);
g.activeAreaVerified=false;
g.note=['Preserved old CAD fill fraction is provisional. Treating the supplied cell footprint as active area ' ...
    'is unverified. A real purchased cell cannot be resized without selecting another physical cell. ' ...
    'Area feasibility alone does not establish rectangular packing, gaps or connector/hinge clearance.'];
end
