function energy = compute_energy_node(node, psm)
    if isempty(node.children)
        energy = [];
        return;
    end

    children = node.children;
    
    for i = 1:length(children)
        children_absx(i) = psm(children(i)).absx;
        children_absy(i) = psm(children(i)).absy;
        children_x(i) = psm(children(i)).x;
        children_y(i) = psm(children(i)).y;
    end
    
    temp = repmat([node.absx; node.absy], 1, length(children)) - [children_absx; children_absy];
    temp = 0.5*sum(temp.^2, 1);
    energy = min(temp);
end