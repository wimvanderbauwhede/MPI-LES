$HOME/Git/RefactorF4Acc/refactorF4acc.pl -P translate_to_OpenCL -c rf4a_to_C.cfg adam_bondv1_feedbf_les_press_v_etc_superkernel 
cp  module_adam_bondv1_feedbf_les_press_v_etc_superkernel.cl module_adam_bondv1_feedbf_les_press_v_etc_superkernel_ORIG.cl
cpp -DNTH=1 -DNUNITS=8 -I. -P module_adam_bondv1_feedbf_les_press_v_etc_superkernel.cl > tmpp.cl
mv tmpp.cl module_adam_bondv1_feedbf_les_press_v_etc_superkernel_after_CPP_for_debugging.cl