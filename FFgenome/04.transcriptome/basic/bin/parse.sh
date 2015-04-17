awk '{if ($4 < 2) print $8"\t"$7"\t"$9"\t"$4}' Root_R8_L1_1.fa.soapout > Root1.soapout
awk '{if ($4 < 2) print $8"\t"$7"\t"$9"\t"$4}' shoot_2_L1_1.fa.soapout > shoot1.soapout
awk '{if ($4 < 2) print $8"\t"$7"\t"$9"\t"$4}' Root_R8_L1_2.fa.soapout > Root2.soapout
awk '{if ($4 < 2) print $8"\t"$7"\t"$9"\t"$4}' shoot_2_L1_2.fa.soapout > shoot2.soapout


awk '{if($4<2) print $8"\t"$7"\t"$1"\t"$1+$6"\t"$4}' Root_R8_L1_1.fa.soapout > root1.soap
