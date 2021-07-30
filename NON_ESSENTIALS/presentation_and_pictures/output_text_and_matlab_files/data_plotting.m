function data_plotting(peptide, score_list, distance_list)
    minREU = min(score_list);
    for i=(1:size(score_list,2))
        if score_list(i) == minREU
            index = i;
        end
    end
    figure
    z = scatter(distance_list,score_list,'h','filled');
    xlim([0 30])
    ylim([-1000 2000])
    z.LineWidth = 0.6;
    z.MarkerEdgeColor = 'k';
    z.MarkerFaceColor = [0 0 1];
    xlabel("Distance [Target = 3.8 A]");
    ylabel("REU (with constraints)");
    title_name = "Enzyme:C2GnT-L, Peptide: " + peptide;
    title(title_name)
    hold on;
    s = scatter(distance_list(index), score_list(index),75,'r','h','filled');
    s.MarkerEdgeColor = 'k';
    output_name = "3OTK_" + peptide + "_funnel";
    saveas(z,output_name,"jpg")
end
