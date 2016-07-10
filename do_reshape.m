function out = do_reshape(row, nrows, ncols, vertices)

    out = [];
    for(kidx = 1: vertices)
        temp_row = reshape(row{kidx}, 1, nrows * ncols);
        
        out = [out temp_row];
    end
end

