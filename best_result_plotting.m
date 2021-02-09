function best_result_plotting(peptide, score_list, distance_list)
    minREU = min(score_list);
    for i=(1:size(score_list,2))
        if score_list(i) == minREU
            index = i;
        end
    end
    scatter(distance_list(index), score_list(index),'MarkerFaceColor',[rand rand rand]);
    xlabel("Distance [Target = 3.8 A]");
    ylabel("REU (with constraints)");
    hold on
end
