echo "OBa TE methylation"
perl drawChrFeature.pl --chrlen ../input/OBa.chrlen --TEfeature ../input/OBa_feature_wave --methylation ../input/OBa_methylation_wave --project OBa > log 2> log2 &

echo "OBa manual TE methylation"
perl drawChrFeature.pl --chrlen ../input/OBa.chrlen --TEfeature ../input/OBa_manual_feature_wave --methylation ../input/OBa_methylation_wave --project OBa > log 2> log2 &

echo "OBa manual TE absolute methylation"
perl drawChrFeature.pl --chrlen ../input/OBa.chrlen --TEfeature ../input/OBa_manual_feature_wave --methylation ../input/OBa_absolute_methylation_wave --project OBa > log 2> log2 &

echo "OBa TE absolute methylation"
perl drawChrFeature.pl --chrlen ../input/OBa.chrlen --TEfeature ../input/OBa_feature_wave --methylation ../input/OBa_absolute_methylation_wave --project OBa > log 2> log2 &

echo "IRGSP TE methylation"
perl drawChrFeature.pl --chrlen ../input/IRGSP.chrlen --TEfeature ../input/IRGSP_feature_wave --methylation ../input/rice_methylation_wave --project IRGSP > log2 > log2 &

echo "Add heterchromatin on chromosome"
perl drawChrFeature.pl --chrlen ../input/OBa.chrlen --TEfeature ../input/OBa_manual_feature_wave --methylation ../input/OBa_methylation_wave --segment ../input/OBa_methylation_wave --project OBa > log 2> log2 &
perl drawChrFeature.pl --chrlen ../input/IRGSP.chrlen --TEfeature ../input/IRGSP_feature_wave --methylation ../input/rice_methylation_wave --segment ../input/rice_methylation_wave --project IRGSP > log2 > log2 &

