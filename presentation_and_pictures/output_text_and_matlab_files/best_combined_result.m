function best_combined_result(core1_rmsd_list, distance_list, peptide_rmsd_list, score_list, peptide)
    minREU = min(score_list);
    for i=(1:size(score_list,2))
        if score_list(i) == minREU
            index = i;
        end
    end
    scatter(core1_rmsd_list(index), score_list(index),'MarkerFaceColor',[rand rand rand]);
    xlabel("Core1 Sugar RMSD");
    ylabel("REU (with constraints)");
    hold on
end