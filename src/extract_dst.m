function DST = extract_dst(DST,irow)
% extract one line (a record for one observation point) from DST

%fprintf(1,'In %s\n',mfilename);

fn=fieldnames(DST);
for i=1:numel(fn)
    F1=getfield(DST,fn{i});
    %            S = SETFIELD(S,{i,j},'field',{k},V) is equivalent to the syntax
    %         S(i,j).field(k) = V;
    %
    %F1(irow)
    
    %fprintf(1,'DST field named %40s %20.10e\n',char(fn{i}),F1(irow));
    DST=setfield(DST,fn{i},F1(irow));
end
return
end

