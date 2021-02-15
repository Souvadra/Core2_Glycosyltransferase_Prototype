function myFunction(core1_rmsd_list, distance_list, peptide_rmsd_list, score_list,peptide) 
    minREU = min(score_list);
    for i=(1:size(score_list,2))
        if score_list(i) == minREU
            index = i;
        end
    end
    figure
    z = scatter(distance_list, score_list, 'h', 'filled');
    xlabel("Distance [Target = 3.8 A]");
    ylabel("REU (with constraints)");
    title_name = "Enzyme:C2GnT-L, Peptide: " + peptide;
    title(title_name)
    hold on;
    s = scatter(distance_list(index), score_list(index),75,'r','h','filled');
    s.MarkerEdgeColor = 'k';
    output_name = "3OTK_" + peptide + "_distance_funnel";
    saveas(z,output_name,"jpg")
    
    figure
    z2 = scatter(core1_rmsd_list, score_list, 'd', 'filled');
    xlabel("RMSD of the Core1 sugar");
    ylabel("REU (with constraints)");
    title(title_name)
    hold on;
    s2 = scatter(core1_rmsd_list(index), score_list(index),75,'r','d','filled');
    s2.MarkerEdgeColor = 'k';
    output_name = "3OTK_" + peptide + "_core1_rmsd_funnel";
    saveas(z2,output_name,"jpg")
    
    figure
    z3 = scatter(peptide_rmsd_list, score_list, 'o', 'filled');
    xlabel("RMSD of the peptide sugar");
    ylabel("REU (with constraints)");
    title(title_name)
    hold on;
    s3 = scatter(peptide_rmsd_list(index), score_list(index),75,'r','o','filled');
    s3.MarkerEdgeColor = 'k';
    output_name = "3OTK_" + peptide + "_peptide_rmsd_funnel";
    saveas(z3,output_name,"jpg")   
end