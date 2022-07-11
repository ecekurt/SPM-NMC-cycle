oneyear=years(1).days(1:150);
%%
nmatrx=[];

for i=1:size(oneyear,2)
nmatrx=[nmatrx;oneyear(i).DChr.CLoss;oneyear(i).Chr.CLoss];
end
%%

Tcell=[];
for i=1:size(oneyear,2)
Tcell=[Tcell; oneyear(i).DChr.Tcell; oneyear(i).Chr.Tcell];
end
% for k=1:size(years(i).days,2)
%     nmatrx=[nmatrx;years(i).days(k).DChr.CLoss;years(i).days(k).Chr.CLoss];
% end