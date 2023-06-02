CAMI= readtable('kallisto_abundance_tpm_filt_stat.csv','FileType', 'text', 'Delimiter', ',');

false_negatives= CAMI{:, [2, 3, 6]};
false_positives=CAMI{:, [2, 3, 7]};

false_negatives(:,3)= false_negatives(:,3)*100;
false_positives(:,3)=false_positives(:,3)*100;
false_negatives(:,2)= false_negatives(:,2)*100;
false_positives(:,2)=false_positives(:,2)*100;

scatter3(false_negatives(:,1), false_negatives(:, 2), false_negatives(:, 3))
xlabel('TPM Threshold')
xlim([0 100])
ylim([0 100])
ylabel('Prevalence (% of samples)')
zlabel('Percent of genomes false positive or false negative')
hold on
scatter3(false_positives(:,1), false_positives(:, 2), false_positives(:, 3), '*')

