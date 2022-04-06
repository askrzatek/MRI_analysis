function newfiles = save_as(namestochange, newnames, destination)

% need unzip files or volumes 
% save files or dir with new names in destination dir if exist or in the
% same folder of namestochange

move = 1;
if ~exist('destination','var')
    move = 0
end


newnames = cellstr(char(newnames));
if ~iscell(namestochange)
    namestochange = cellstr(char(namestochange));
end 
 
if length(newnames) == 1
    mynewnames = cell(length(namestochange),1);
    [mynewnames{:}] =deal(newnames{1});
elseif length(namestochange) ~= length(newnames) & length(newnames) ~= 1
    error('la longeur des cells n''est pas la même');
else
    mynewnames = newnames;
end

[pathnames, allnamestochange] = get_parent_path(namestochange,1);
[~,allnewnames] = get_parent_path(mynewnames,1);

for nbr = 1 :length(allnamestochange)
    extension = '';
    if isfile(namestochange{nbr})
        [~,~,extension] = fileparts( deblank( namestochange{nbr}));
    end
    
    name = fullfile(pathnames{nbr},[allnewnames{nbr} extension]);
    newlist{nbr} = name;

end

listnewnames =newlist';



if move
    r_movefile(namestochange,listnewnames,'copy');
    newfiles = r_movefile(listnewnames,destination,'move');
else
    
    newfiles  = r_movefile(namestochange,listnewnames,'copy');
end
end