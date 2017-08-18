$HOME/Git/RefactorF4Acc/refactorF4acc.pl -P translate_to_OpenCL -c rf4a_to_C.cfg adam_bondv1_feedbf_les_press_v_etc_superkernel 
cp  module_adam_bondv1_feedbf_les_press_v_etc_superkernel.cl module_adam_bondv1_feedbf_les_press_v_etc_superkernel_ORIG.cl
MACROS_KERNEL=`perl -e '@ls=<>;$str=join(" ",map {$_=~s/\n//;s/^\s*//;s/\s*$//;s/.define\s*/-D/;s/.undef\s*/-U/;s/\s+/=/;$_} @ls);print $str' ./macros_kernel.h`
echo $MACROS_KERNEL
cpp $MACROS_KERNEL -I. -P module_adam_bondv1_feedbf_les_press_v_etc_superkernel.cl > module_adam_bondv1_feedbf_les_press_v_etc_superkernel_after_CPP_for_debugging.cl
