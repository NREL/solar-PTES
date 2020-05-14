function [] = print_hexs(HX,iL)
%PRINT_HEXS Print a table summarising the main heat exchanger parameters

rows   = length(HX);

name   = {HX.name}';
fHname = cell(rows,1);
fCname = cell(rows,1);
DppH   = zeros(rows,1);
DppC   = zeros(rows,1);
DTmin  = zeros(rows,1);
effDT  = zeros(rows,1);
L      = [HX.L]';
A      = [HX.A1]';
Af     = [HX.Af1]';
COST   = zeros(rows,1);
Sirr   = zeros(rows,1);
for i = 1:length(HX)
    fHname{i} = valid_name(HX(i).H(iL).name,2);
    fCname{i} = valid_name(HX(i).C(iL).name,2);
    DppH(i)   = HX(i).DppH(iL);
    DppC(i)   = HX(i).DppC(iL);
    DTmin(i)  = HX(i).DTmin(iL);
    effDT(i)  = HX(i).effDT(iL);
    COST(i)   = HX(i).hx_cost.COST;
    Sirr(i)   = sum(HX(i).Sirr(iL,:));
end

%HX_summary = table(name,fHname,fCname,DTmin,effDT,DppH,DppC,L,A,Af);
%disp(HX_summary)

fprintf('%10s ','Name','fHname','fCname','DTmin','effDT','DppH','DppC','L','A','Af','COST','Sirr')
fprintf('\n')
for i = 1:length(HX)
    fprintf('%10s %10s %10s %10.2f %10.3f %10.3f %10.3f %10.1f %10.2e %10.1f %10.1e %10.1e\n',...
        name{i},fHname{i},fCname{i},DTmin(i),effDT(i),DppH(i),DppC(i),L(i),A(i),Af(i),COST(i),Sirr(i))
end
fprintf('\n')

end

