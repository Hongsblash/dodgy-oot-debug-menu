glabel func_80AE067C
/* 0292C 80AE067C 8C820190 */  lw      $v0, 0x0190($a0)           ## 00000190
/* 02930 80AE0680 3C0E80AE */  lui     $t6, %hi(func_80ADEDFC)    ## $t6 = 80AE0000
/* 02934 80AE0684 25CEEDFC */  addiu   $t6, $t6, %lo(func_80ADEDFC) ## $t6 = 80ADEDFC
/* 02938 80AE0688 15C20021 */  bne     $t6, $v0, .L80AE0710       
/* 0293C 80AE068C 3C1980AE */  lui     $t9, %hi(func_80ADF894)    ## $t9 = 80AE0000
/* 02940 80AE0690 9082029A */  lbu     $v0, 0x029A($a0)           ## 0000029A
/* 02944 80AE0694 240F00FF */  addiu   $t7, $zero, 0x00FF         ## $t7 = 000000FF
/* 02948 80AE0698 24180032 */  addiu   $t8, $zero, 0x0032         ## $t8 = 00000032
/* 0294C 80AE069C 24420005 */  addiu   $v0, $v0, 0x0005           ## $v0 = 00000005
/* 02950 80AE06A0 00021400 */  sll     $v0, $v0, 16               
/* 02954 80AE06A4 00021403 */  sra     $v0, $v0, 16               
/* 02958 80AE06A8 28410100 */  slti    $at, $v0, 0x0100           
/* 0295C 80AE06AC 54200004 */  bnel    $at, $zero, .L80AE06C0     
/* 02960 80AE06B0 A082029A */  sb      $v0, 0x029A($a0)           ## 0000029A
/* 02964 80AE06B4 10000002 */  beq     $zero, $zero, .L80AE06C0   
/* 02968 80AE06B8 A08F029A */  sb      $t7, 0x029A($a0)           ## 0000029A
/* 0296C 80AE06BC A082029A */  sb      $v0, 0x029A($a0)           ## 0000029A
.L80AE06C0:
/* 02970 80AE06C0 9082029B */  lbu     $v0, 0x029B($a0)           ## 0000029B
/* 02974 80AE06C4 2442FFFB */  addiu   $v0, $v0, 0xFFFB           ## $v0 = 00000000
/* 02978 80AE06C8 00021400 */  sll     $v0, $v0, 16               
/* 0297C 80AE06CC 00021403 */  sra     $v0, $v0, 16               
/* 02980 80AE06D0 28410032 */  slti    $at, $v0, 0x0032           
/* 02984 80AE06D4 50200004 */  beql    $at, $zero, .L80AE06E8     
/* 02988 80AE06D8 A082029B */  sb      $v0, 0x029B($a0)           ## 0000029B
/* 0298C 80AE06DC 10000002 */  beq     $zero, $zero, .L80AE06E8   
/* 02990 80AE06E0 A098029B */  sb      $t8, 0x029B($a0)           ## 0000029B
/* 02994 80AE06E4 A082029B */  sb      $v0, 0x029B($a0)           ## 0000029B
.L80AE06E8:
/* 02998 80AE06E8 9082029C */  lbu     $v0, 0x029C($a0)           ## 0000029C
/* 0299C 80AE06EC 2442FFFB */  addiu   $v0, $v0, 0xFFFB           ## $v0 = FFFFFFFB
/* 029A0 80AE06F0 00021400 */  sll     $v0, $v0, 16               
/* 029A4 80AE06F4 00021403 */  sra     $v0, $v0, 16               
/* 029A8 80AE06F8 04410003 */  bgez    $v0, .L80AE0708            
/* 029AC 80AE06FC 00000000 */  nop
/* 029B0 80AE0700 03E00008 */  jr      $ra                        
/* 029B4 80AE0704 A080029C */  sb      $zero, 0x029C($a0)         ## 0000029C
.L80AE0708:
/* 029B8 80AE0708 03E00008 */  jr      $ra                        
/* 029BC 80AE070C A082029C */  sb      $v0, 0x029C($a0)           ## 0000029C
.L80AE0710:
/* 029C0 80AE0710 2739F894 */  addiu   $t9, $t9, %lo(func_80ADF894) ## $t9 = 80ADF894
/* 029C4 80AE0714 17220022 */  bne     $t9, $v0, .L80AE07A0       
/* 029C8 80AE0718 3C0B80AE */  lui     $t3, %hi(func_80ADEECC)    ## $t3 = 80AE0000
/* 029CC 80AE071C 9082029A */  lbu     $v0, 0x029A($a0)           ## 0000029A
/* 029D0 80AE0720 24080050 */  addiu   $t0, $zero, 0x0050         ## $t0 = 00000050
/* 029D4 80AE0724 240900FF */  addiu   $t1, $zero, 0x00FF         ## $t1 = 000000FF
/* 029D8 80AE0728 24420005 */  addiu   $v0, $v0, 0x0005           ## $v0 = 00000000
/* 029DC 80AE072C 00021400 */  sll     $v0, $v0, 16               
/* 029E0 80AE0730 00021403 */  sra     $v0, $v0, 16               
/* 029E4 80AE0734 28410051 */  slti    $at, $v0, 0x0051           
/* 029E8 80AE0738 14200003 */  bne     $at, $zero, .L80AE0748     
/* 029EC 80AE073C 240A00E1 */  addiu   $t2, $zero, 0x00E1         ## $t2 = 000000E1
/* 029F0 80AE0740 10000002 */  beq     $zero, $zero, .L80AE074C   
/* 029F4 80AE0744 A088029A */  sb      $t0, 0x029A($a0)           ## 0000029A
.L80AE0748:
/* 029F8 80AE0748 A082029A */  sb      $v0, 0x029A($a0)           ## 0000029A
.L80AE074C:
/* 029FC 80AE074C 9082029B */  lbu     $v0, 0x029B($a0)           ## 0000029B
/* 02A00 80AE0750 24420005 */  addiu   $v0, $v0, 0x0005           ## $v0 = 00000005
/* 02A04 80AE0754 00021400 */  sll     $v0, $v0, 16               
/* 02A08 80AE0758 00021403 */  sra     $v0, $v0, 16               
/* 02A0C 80AE075C 28410100 */  slti    $at, $v0, 0x0100           
/* 02A10 80AE0760 54200004 */  bnel    $at, $zero, .L80AE0774     
/* 02A14 80AE0764 A082029B */  sb      $v0, 0x029B($a0)           ## 0000029B
/* 02A18 80AE0768 10000002 */  beq     $zero, $zero, .L80AE0774   
/* 02A1C 80AE076C A089029B */  sb      $t1, 0x029B($a0)           ## 0000029B
/* 02A20 80AE0770 A082029B */  sb      $v0, 0x029B($a0)           ## 0000029B
.L80AE0774:
/* 02A24 80AE0774 9082029C */  lbu     $v0, 0x029C($a0)           ## 0000029C
/* 02A28 80AE0778 24420005 */  addiu   $v0, $v0, 0x0005           ## $v0 = 0000000A
/* 02A2C 80AE077C 00021400 */  sll     $v0, $v0, 16               
/* 02A30 80AE0780 00021403 */  sra     $v0, $v0, 16               
/* 02A34 80AE0784 284100E2 */  slti    $at, $v0, 0x00E2           
/* 02A38 80AE0788 14200003 */  bne     $at, $zero, .L80AE0798     
/* 02A3C 80AE078C 00000000 */  nop
/* 02A40 80AE0790 03E00008 */  jr      $ra                        
/* 02A44 80AE0794 A08A029C */  sb      $t2, 0x029C($a0)           ## 0000029C
.L80AE0798:
/* 02A48 80AE0798 03E00008 */  jr      $ra                        
/* 02A4C 80AE079C A082029C */  sb      $v0, 0x029C($a0)           ## 0000029C
.L80AE07A0:
/* 02A50 80AE07A0 256BEECC */  addiu   $t3, $t3, %lo(func_80ADEECC) ## $t3 = 80ADEECC
/* 02A54 80AE07A4 55620010 */  bnel    $t3, $v0, .L80AE07E8       
/* 02A58 80AE07A8 9082029A */  lbu     $v0, 0x029A($a0)           ## 0000029A
/* 02A5C 80AE07AC 908C0114 */  lbu     $t4, 0x0114($a0)           ## 00000114
/* 02A60 80AE07B0 240E0050 */  addiu   $t6, $zero, 0x0050         ## $t6 = 00000050
/* 02A64 80AE07B4 240F00FF */  addiu   $t7, $zero, 0x00FF         ## $t7 = 000000FF
/* 02A68 80AE07B8 318D0002 */  andi    $t5, $t4, 0x0002           ## $t5 = 00000000
/* 02A6C 80AE07BC 11A00005 */  beq     $t5, $zero, .L80AE07D4     
/* 02A70 80AE07C0 241800E1 */  addiu   $t8, $zero, 0x00E1         ## $t8 = 000000E1
/* 02A74 80AE07C4 A080029A */  sb      $zero, 0x029A($a0)         ## 0000029A
/* 02A78 80AE07C8 A080029B */  sb      $zero, 0x029B($a0)         ## 0000029B
/* 02A7C 80AE07CC 03E00008 */  jr      $ra                        
/* 02A80 80AE07D0 A080029C */  sb      $zero, 0x029C($a0)         ## 0000029C
.L80AE07D4:
/* 02A84 80AE07D4 A08E029A */  sb      $t6, 0x029A($a0)           ## 0000029A
/* 02A88 80AE07D8 A08F029B */  sb      $t7, 0x029B($a0)           ## 0000029B
/* 02A8C 80AE07DC 03E00008 */  jr      $ra                        
/* 02A90 80AE07E0 A098029C */  sb      $t8, 0x029C($a0)           ## 0000029C
.L80AE07E4:
/* 02A94 80AE07E4 9082029A */  lbu     $v0, 0x029A($a0)           ## 0000029A
.L80AE07E8:
/* 02A98 80AE07E8 241900FF */  addiu   $t9, $zero, 0x00FF         ## $t9 = 000000FF
/* 02A9C 80AE07EC 240800FF */  addiu   $t0, $zero, 0x00FF         ## $t0 = 000000FF
/* 02AA0 80AE07F0 24420005 */  addiu   $v0, $v0, 0x0005           ## $v0 = 0000000F
/* 02AA4 80AE07F4 00021400 */  sll     $v0, $v0, 16               
/* 02AA8 80AE07F8 00021403 */  sra     $v0, $v0, 16               
/* 02AAC 80AE07FC 28410100 */  slti    $at, $v0, 0x0100           
/* 02AB0 80AE0800 54200004 */  bnel    $at, $zero, .L80AE0814     
/* 02AB4 80AE0804 A082029A */  sb      $v0, 0x029A($a0)           ## 0000029A
/* 02AB8 80AE0808 10000002 */  beq     $zero, $zero, .L80AE0814   
/* 02ABC 80AE080C A099029A */  sb      $t9, 0x029A($a0)           ## 0000029A
/* 02AC0 80AE0810 A082029A */  sb      $v0, 0x029A($a0)           ## 0000029A
.L80AE0814:
/* 02AC4 80AE0814 9082029B */  lbu     $v0, 0x029B($a0)           ## 0000029B
/* 02AC8 80AE0818 24420005 */  addiu   $v0, $v0, 0x0005           ## $v0 = 00000014
/* 02ACC 80AE081C 00021400 */  sll     $v0, $v0, 16               
/* 02AD0 80AE0820 00021403 */  sra     $v0, $v0, 16               
/* 02AD4 80AE0824 28410100 */  slti    $at, $v0, 0x0100           
/* 02AD8 80AE0828 54200004 */  bnel    $at, $zero, .L80AE083C     
/* 02ADC 80AE082C A082029B */  sb      $v0, 0x029B($a0)           ## 0000029B
/* 02AE0 80AE0830 10000002 */  beq     $zero, $zero, .L80AE083C   
/* 02AE4 80AE0834 A088029B */  sb      $t0, 0x029B($a0)           ## 0000029B
/* 02AE8 80AE0838 A082029B */  sb      $v0, 0x029B($a0)           ## 0000029B
.L80AE083C:
/* 02AEC 80AE083C 9083029C */  lbu     $v1, 0x029C($a0)           ## 0000029C
/* 02AF0 80AE0840 286100D3 */  slti    $at, $v1, 0x00D3           
/* 02AF4 80AE0844 1420000B */  bne     $at, $zero, .L80AE0874     
/* 02AF8 80AE0848 24620005 */  addiu   $v0, $v1, 0x0005           ## $v0 = 00000005
/* 02AFC 80AE084C 2462FFFB */  addiu   $v0, $v1, 0xFFFB           ## $v0 = FFFFFFFB
/* 02B00 80AE0850 00021400 */  sll     $v0, $v0, 16               
/* 02B04 80AE0854 00021403 */  sra     $v0, $v0, 16               
/* 02B08 80AE0858 284100D2 */  slti    $at, $v0, 0x00D2           
/* 02B0C 80AE085C 10200003 */  beq     $at, $zero, .L80AE086C     
/* 02B10 80AE0860 240900D2 */  addiu   $t1, $zero, 0x00D2         ## $t1 = 000000D2
/* 02B14 80AE0864 03E00008 */  jr      $ra                        
/* 02B18 80AE0868 A089029C */  sb      $t1, 0x029C($a0)           ## 0000029C
.L80AE086C:
/* 02B1C 80AE086C 03E00008 */  jr      $ra                        
/* 02B20 80AE0870 A082029C */  sb      $v0, 0x029C($a0)           ## 0000029C
.L80AE0874:
/* 02B24 80AE0874 00021400 */  sll     $v0, $v0, 16               
/* 02B28 80AE0878 00021403 */  sra     $v0, $v0, 16               
/* 02B2C 80AE087C 284100D3 */  slti    $at, $v0, 0x00D3           
/* 02B30 80AE0880 14200003 */  bne     $at, $zero, .L80AE0890     
/* 02B34 80AE0884 240A00D2 */  addiu   $t2, $zero, 0x00D2         ## $t2 = 000000D2
/* 02B38 80AE0888 03E00008 */  jr      $ra                        
/* 02B3C 80AE088C A08A029C */  sb      $t2, 0x029C($a0)           ## 0000029C
.L80AE0890:
/* 02B40 80AE0890 A082029C */  sb      $v0, 0x029C($a0)           ## 0000029C
/* 02B44 80AE0894 03E00008 */  jr      $ra                        
/* 02B48 80AE0898 00000000 */  nop


