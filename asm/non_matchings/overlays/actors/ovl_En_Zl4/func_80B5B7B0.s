glabel func_80B5B7B0
/* 00000 80B5B7B0 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 00004 80B5B7B4 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00008 80B5B7B8 AFA5002C */  sw      $a1, 0x002C($sp)           
/* 0000C 80B5B7BC 848E07A0 */  lh      $t6, 0x07A0($a0)           ## 000007A0
/* 00010 80B5B7C0 24050021 */  addiu   $a1, $zero, 0x0021         ## $a1 = 00000021
/* 00014 80B5B7C4 000E7880 */  sll     $t7, $t6,  2               
/* 00018 80B5B7C8 008FC021 */  addu    $t8, $a0, $t7              
/* 0001C 80B5B7CC 8F040790 */  lw      $a0, 0x0790($t8)           ## 00000790
/* 00020 80B5B7D0 0C0169DF */  jal     func_8005A77C              
/* 00024 80B5B7D4 AFA40024 */  sw      $a0, 0x0024($sp)           
/* 00028 80B5B7D8 87B9002E */  lh      $t9, 0x002E($sp)           
/* 0002C 80B5B7DC 3C0980B6 */  lui     $t1, %hi(D_80B5EAE8)       ## $t1 = 80B60000
/* 00030 80B5B7E0 2529EAE8 */  addiu   $t1, $t1, %lo(D_80B5EAE8)  ## $t1 = 80B5EAE8
/* 00034 80B5B7E4 001940C0 */  sll     $t0, $t9,  3               
/* 00038 80B5B7E8 01194023 */  subu    $t0, $t0, $t9              
/* 0003C 80B5B7EC 00084080 */  sll     $t0, $t0,  2               
/* 00040 80B5B7F0 01091021 */  addu    $v0, $t0, $t1              
/* 00044 80B5B7F4 8FA40024 */  lw      $a0, 0x0024($sp)           
/* 00048 80B5B7F8 8C4B0000 */  lw      $t3, 0x0000($v0)           ## 00000000
/* 0004C 80B5B7FC 27A30018 */  addiu   $v1, $sp, 0x0018           ## $v1 = FFFFFFF0
/* 00050 80B5B800 AC8B0050 */  sw      $t3, 0x0050($a0)           ## 00000050
/* 00054 80B5B804 8C4A0004 */  lw      $t2, 0x0004($v0)           ## 00000004
/* 00058 80B5B808 AC8A0054 */  sw      $t2, 0x0054($a0)           ## 00000054
/* 0005C 80B5B80C 8C4B0008 */  lw      $t3, 0x0008($v0)           ## 00000008
/* 00060 80B5B810 AC8B0058 */  sw      $t3, 0x0058($a0)           ## 00000058
/* 00064 80B5B814 8C4D000C */  lw      $t5, 0x000C($v0)           ## 0000000C
/* 00068 80B5B818 8C4C0010 */  lw      $t4, 0x0010($v0)           ## 00000010
/* 0006C 80B5B81C AC6D0000 */  sw      $t5, 0x0000($v1)           ## FFFFFFF0
/* 00070 80B5B820 8C4D0014 */  lw      $t5, 0x0014($v0)           ## 00000014
/* 00074 80B5B824 8C6F0000 */  lw      $t7, 0x0000($v1)           ## FFFFFFF0
/* 00078 80B5B828 AC6C0004 */  sw      $t4, 0x0004($v1)           ## FFFFFFF4
/* 0007C 80B5B82C AC6D0008 */  sw      $t5, 0x0008($v1)           ## FFFFFFF8
/* 00080 80B5B830 AC8F0074 */  sw      $t7, 0x0074($a0)           ## 00000074
/* 00084 80B5B834 8C6E0004 */  lw      $t6, 0x0004($v1)           ## FFFFFFF4
/* 00088 80B5B838 AC8E0078 */  sw      $t6, 0x0078($a0)           ## 00000078
/* 0008C 80B5B83C 8C6F0008 */  lw      $t7, 0x0008($v1)           ## FFFFFFF8
/* 00090 80B5B840 AC8F007C */  sw      $t7, 0x007C($a0)           ## 0000007C
/* 00094 80B5B844 8C790000 */  lw      $t9, 0x0000($v1)           ## FFFFFFF0
/* 00098 80B5B848 AC99005C */  sw      $t9, 0x005C($a0)           ## 0000005C
/* 0009C 80B5B84C 8C780004 */  lw      $t8, 0x0004($v1)           ## FFFFFFF4
/* 000A0 80B5B850 AC980060 */  sw      $t8, 0x0060($a0)           ## 00000060
/* 000A4 80B5B854 8C790008 */  lw      $t9, 0x0008($v1)           ## FFFFFFF8
/* 000A8 80B5B858 AC990064 */  sw      $t9, 0x0064($a0)           ## 00000064
/* 000AC 80B5B85C 84480018 */  lh      $t0, 0x0018($v0)           ## 00000018
/* 000B0 80B5B860 A488015A */  sh      $t0, 0x015A($a0)           ## 0000015A
/* 000B4 80B5B864 8449001A */  lh      $t1, 0x001A($v0)           ## 0000001A
/* 000B8 80B5B868 44892000 */  mtc1    $t1, $f4                   ## $f4 = -0.00
/* 000BC 80B5B86C 00000000 */  nop
/* 000C0 80B5B870 468021A0 */  cvt.s.w $f6, $f4                   
/* 000C4 80B5B874 E48600FC */  swc1    $f6, 0x00FC($a0)           ## 000000FC
/* 000C8 80B5B878 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 000CC 80B5B87C 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 000D0 80B5B880 03E00008 */  jr      $ra                        
/* 000D4 80B5B884 00000000 */  nop


