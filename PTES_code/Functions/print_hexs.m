function [] = print_hexs(HX,iL,varargin)
%PRINT_HEXS Print a table summarising the main heat exchanger parameters

if isempty(iL)
    return
else
    iL = iL(1);
end

% Optionally, print text contained in varargin (e.g. 'Charge cycle:')
switch nargin
    case 2
    case 3
        fprintf([varargin{1},'\n'])
    otherwise
        error('not implemented')
end

% Preallocate arrays of fields to be printed
rows   = length(HX);
name   = cell(rows,1);
fHname = cell(rows,1);
fCname = cell(rows,1);
DppH   = zeros(rows,1);
DppC   = zeros(rows,1);
DTmin  = zeros(rows,1);
effDT  = zeros(rows,1);
L1      = zeros(rows,1);
A1      = zeros(rows,1);
Af     = zeros(rows,1);
COST   = zeros(rows,1);
Sirr   = zeros(rows,1);

% Gather content of arrays
for i = 1:length(HX)
    name{i}   = HX(i).name;
    fHname{i} = valid_name(HX(i).H(iL).name,2);
    fCname{i} = valid_name(HX(i).C(iL).name,2);
    DppH(i)   = HX(i).DppH(iL);
    DppC(i)   = HX(i).DppC(iL);
    DTmin(i)  = HX(i).DTmin(iL);
    effDT(i)  = HX(i).effDT(iL);
    if~isempty(HX(i).L1)
        L1(i) = HX(i).L1;
        A1(i) = HX(i).A1;
        Af(i) = HX(i).Af1;
    end
    if~isempty(HX(i).hx_cost.COST)
        COST(i) = HX(i).hx_cost.COST;
    end
    Sirr(i)   = sum(HX(i).Sirr(iL,:));
end


% Print header
fprintf('%10s ','Name','fHname','fCname','DTmin','effDT','DppH','DppC','L','A','Af','COST','Sirr')
fprintf('\n')

% Print a summary of each heat exchanger, one by one
for i = 1:length(HX)
    if ~isempty(fHname{i})
        fprintf('%10s %10s %10s %10.2f %10.3f %10.3f %10.3f %10.1f %10.2e %10.1f %10.1e %10.1e\n',...
            name{i},fHname{i},fCname{i},DTmin(i),effDT(i),DppH(i),DppC(i),L1(i),A1(i),Af(i),COST(i),Sirr(i))
    end
end
fprintf('\n')

% Alternative way of printing list of HXs using Matlab's table structure
%HX_summary = table(name,fHname,fCname,DTmin,effDT,DppH,DppC,L,A,Af,COST,Sirr);
%disp(HX_summary)

end

