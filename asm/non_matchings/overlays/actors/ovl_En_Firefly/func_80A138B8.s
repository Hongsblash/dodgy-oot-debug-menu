glabel func_80A138B8
/* 00848 80A138B8 27BDFFB0 */  addiu   $sp, $sp, 0xFFB0           ## $sp = FFFFFFB0
/* 0084C 80A138BC AFBF002C */  sw      $ra, 0x002C($sp)           
/* 00850 80A138C0 AFB30028 */  sw      $s3, 0x0028($sp)           
/* 00854 80A138C4 AFB20024 */  sw      $s2, 0x0024($sp)           
/* 00858 80A138C8 AFB10020 */  sw      $s1, 0x0020($sp)           
/* 0085C 80A138CC AFB0001C */  sw      $s0, 0x001C($sp)           
/* 00860 80A138D0 F7B40010 */  sdc1    $f20, 0x0010($sp)          
/* 00864 80A138D4 8CB01C64 */  lw      $s0, 0x1C64($a1)           ## 00001C64
/* 00868 80A138D8 3C0180A1 */  lui     $at, %hi(D_80A1503C)       ## $at = 80A10000
/* 0086C 80A138DC 00809825 */  or      $s3, $a0, $zero            ## $s3 = 00000000
/* 00870 80A138E0 00008825 */  or      $s1, $zero, $zero          ## $s1 = 00000000
/* 00874 80A138E4 12000014 */  beq     $s0, $zero, .L80A13938     
/* 00878 80A138E8 C434503C */  lwc1    $f20, %lo(D_80A1503C)($at) 
/* 0087C 80A138EC 2412005E */  addiu   $s2, $zero, 0x005E         ## $s2 = 0000005E
/* 00880 80A138F0 860E0000 */  lh      $t6, 0x0000($s0)           ## 00000000
.L80A138F4:
/* 00884 80A138F4 564E000E */  bnel    $s2, $t6, .L80A13930       
/* 00888 80A138F8 8E100124 */  lw      $s0, 0x0124($s0)           ## 00000124
/* 0088C 80A138FC 860F01E4 */  lh      $t7, 0x01E4($s0)           ## 000001E4
/* 00890 80A13900 02602025 */  or      $a0, $s3, $zero            ## $a0 = 00000000
/* 00894 80A13904 51E0000A */  beql    $t7, $zero, .L80A13930     
/* 00898 80A13908 8E100124 */  lw      $s0, 0x0124($s0)           ## 00000124
/* 0089C 80A1390C 0C00B6D2 */  jal     func_8002DB48              
/* 008A0 80A13910 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 008A4 80A13914 4614003C */  c.lt.s  $f0, $f20                  
/* 008A8 80A13918 00000000 */  nop
/* 008AC 80A1391C 45020004 */  bc1fl   .L80A13930                 
/* 008B0 80A13920 8E100124 */  lw      $s0, 0x0124($s0)           ## 00000124
/* 008B4 80A13924 46000506 */  mov.s   $f20, $f0                  
/* 008B8 80A13928 02008825 */  or      $s1, $s0, $zero            ## $s1 = 00000000
/* 008BC 80A1392C 8E100124 */  lw      $s0, 0x0124($s0)           ## 00000124
.L80A13930:
/* 008C0 80A13930 5600FFF0 */  bnel    $s0, $zero, .L80A138F4     
/* 008C4 80A13934 860E0000 */  lh      $t6, 0x0000($s0)           ## 00000000
.L80A13938:
/* 008C8 80A13938 1220002B */  beq     $s1, $zero, .L80A139E8     
/* 008CC 80A1393C 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 008D0 80A13940 C6240024 */  lwc1    $f4, 0x0024($s1)           ## 00000024
/* 008D4 80A13944 3C014170 */  lui     $at, 0x4170                ## $at = 41700000
/* 008D8 80A13948 4481A000 */  mtc1    $at, $f20                  ## $f20 = 15.00
/* 008DC 80A1394C 3C014250 */  lui     $at, 0x4250                ## $at = 42500000
/* 008E0 80A13950 E7A40034 */  swc1    $f4, 0x0034($sp)           
/* 008E4 80A13954 C6260028 */  lwc1    $f6, 0x0028($s1)           ## 00000028
/* 008E8 80A13958 44814000 */  mtc1    $at, $f8                   ## $f8 = 52.00
/* 008EC 80A1395C 27B00034 */  addiu   $s0, $sp, 0x0034           ## $s0 = FFFFFFE4
/* 008F0 80A13960 02002825 */  or      $a1, $s0, $zero            ## $a1 = FFFFFFE4
/* 008F4 80A13964 46083280 */  add.s   $f10, $f6, $f8             
/* 008F8 80A13968 02602025 */  or      $a0, $s3, $zero            ## $a0 = 00000000
/* 008FC 80A1396C 46145400 */  add.s   $f16, $f10, $f20           
/* 00900 80A13970 E7B00038 */  swc1    $f16, 0x0038($sp)          
/* 00904 80A13974 C632002C */  lwc1    $f18, 0x002C($s1)          ## 0000002C
/* 00908 80A13978 0C00B6DB */  jal     func_8002DB6C              
/* 0090C 80A1397C E7B2003C */  swc1    $f18, 0x003C($sp)          
/* 00910 80A13980 4614003C */  c.lt.s  $f0, $f20                  
/* 00914 80A13984 02602025 */  or      $a0, $s3, $zero            ## $a0 = 00000000
/* 00918 80A13988 45000005 */  bc1f    .L80A139A0                 
/* 0091C 80A1398C 00000000 */  nop
/* 00920 80A13990 0C284C26 */  jal     func_80A13098              
/* 00924 80A13994 02602025 */  or      $a0, $s3, $zero            ## $a0 = 00000000
/* 00928 80A13998 10000013 */  beq     $zero, $zero, .L80A139E8   
/* 0092C 80A1399C 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L80A139A0:
/* 00930 80A139A0 0C00B69E */  jal     func_8002DA78              
/* 00934 80A139A4 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 00938 80A139A8 00022C00 */  sll     $a1, $v0, 16               
/* 0093C 80A139AC 00052C03 */  sra     $a1, $a1, 16               
/* 00940 80A139B0 266400B6 */  addiu   $a0, $s3, 0x00B6           ## $a0 = 000000B6
/* 00944 80A139B4 0C01DE2B */  jal     Math_ApproxUpdateScaledS
              
/* 00948 80A139B8 24060300 */  addiu   $a2, $zero, 0x0300         ## $a2 = 00000300
/* 0094C 80A139BC 02602025 */  or      $a0, $s3, $zero            ## $a0 = 00000000
/* 00950 80A139C0 0C00B6CA */  jal     func_8002DB28              
/* 00954 80A139C4 02002825 */  or      $a1, $s0, $zero            ## $a1 = FFFFFFE4
/* 00958 80A139C8 24451554 */  addiu   $a1, $v0, 0x1554           ## $a1 = 00001554
/* 0095C 80A139CC 00052C00 */  sll     $a1, $a1, 16               
/* 00960 80A139D0 00052C03 */  sra     $a1, $a1, 16               
/* 00964 80A139D4 266400B4 */  addiu   $a0, $s3, 0x00B4           ## $a0 = 000000B4
/* 00968 80A139D8 0C01DE2B */  jal     Math_ApproxUpdateScaledS
              
/* 0096C 80A139DC 24060100 */  addiu   $a2, $zero, 0x0100         ## $a2 = 00000100
/* 00970 80A139E0 10000001 */  beq     $zero, $zero, .L80A139E8   
/* 00974 80A139E4 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
.L80A139E8:
/* 00978 80A139E8 8FBF002C */  lw      $ra, 0x002C($sp)           
/* 0097C 80A139EC D7B40010 */  ldc1    $f20, 0x0010($sp)          
/* 00980 80A139F0 8FB0001C */  lw      $s0, 0x001C($sp)           
/* 00984 80A139F4 8FB10020 */  lw      $s1, 0x0020($sp)           
/* 00988 80A139F8 8FB20024 */  lw      $s2, 0x0024($sp)           
/* 0098C 80A139FC 8FB30028 */  lw      $s3, 0x0028($sp)           
/* 00990 80A13A00 03E00008 */  jr      $ra                        
/* 00994 80A13A04 27BD0050 */  addiu   $sp, $sp, 0x0050           ## $sp = 00000000


