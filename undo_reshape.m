function out = undo_reshape(row, nrows, ncols, vertices)
    l = length(row) / vertices;
    out = cell(vertices,1);
    
    for (i = 0:(vertices - 1)) 
        start = ((l * i) + 1);
%             param.K{i+1} = ;
        out{i+1} = reshape(row(start:(start+(l - 1))), nrows, ncols);
    end
end

