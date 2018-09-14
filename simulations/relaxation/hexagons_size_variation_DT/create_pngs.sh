for ((B=-2000; B<1; B+=100))
do
    python ~/scripts/mayavi_save_files.py \
    "vtks/2Dhex_hexagon_R6nm_PdFe-Ir_1-sk-up_B${B}mT_vtks" \
    --fidimag --camera_distance 67 \
    --output_name "R20nm_B${B}"
done
