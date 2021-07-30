function All_Data_plotting(peptide)
    run(peptide)
    minREU = min(score_list);
    for i=(1:size(score_list,2))
        if score_list(i) == minREU
            index = i;
        end
    end
    figure
    scatter(distance_list, score_list, 'h', 'filled');
    xlabel("Distance [Target = 3.8 A]");
    ylabel("REU (with constriants)"); 
    title_name = "Enzyme: C2GnT-L, Peptide: " + peptide;
    
end