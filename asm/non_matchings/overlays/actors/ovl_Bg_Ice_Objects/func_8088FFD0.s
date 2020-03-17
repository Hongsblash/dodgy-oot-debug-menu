glabel func_8088FFD0
/* 007C0 8088FFD0 27BDFF98 */  addiu   $sp, $sp, 0xFF98           ## $sp = FFFFFF98
/* 007C4 8088FFD4 AFB00030 */  sw      $s0, 0x0030($sp)           
/* 007C8 8088FFD8 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 007CC 8088FFDC AFBF0034 */  sw      $ra, 0x0034($sp)           
/* 007D0 8088FFE0 AFA5006C */  sw      $a1, 0x006C($sp)           
/* 007D4 8088FFE4 F7B40028 */  sdc1    $f20, 0x0028($sp)          
/* 007D8 8088FFE8 3C054120 */  lui     $a1, 0x4120                ## $a1 = 41200000
/* 007DC 8088FFEC 24840068 */  addiu   $a0, $a0, 0x0068           ## $a0 = 00000068
/* 007E0 8088FFF0 0C01DE80 */  jal     Math_ApproxF
              
/* 007E4 8088FFF4 3C063F00 */  lui     $a2, 0x3F00                ## $a2 = 3F000000
/* 007E8 8088FFF8 26040024 */  addiu   $a0, $s0, 0x0024           ## $a0 = 00000024
/* 007EC 8088FFFC 8E050168 */  lw      $a1, 0x0168($s0)           ## 00000168
/* 007F0 80890000 0C01DE80 */  jal     Math_ApproxF
              
/* 007F4 80890004 8E060068 */  lw      $a2, 0x0068($s0)           ## 00000068
/* 007F8 80890008 8E050170 */  lw      $a1, 0x0170($s0)           ## 00000170
/* 007FC 8089000C 8E060068 */  lw      $a2, 0x0068($s0)           ## 00000068
/* 00800 80890010 AFA20064 */  sw      $v0, 0x0064($sp)           
/* 00804 80890014 0C01DE80 */  jal     Math_ApproxF
              
/* 00808 80890018 2604002C */  addiu   $a0, $s0, 0x002C           ## $a0 = 0000002C
/* 0080C 8089001C 8FA30064 */  lw      $v1, 0x0064($sp)           
/* 00810 80890020 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 00814 80890024 24060007 */  addiu   $a2, $zero, 0x0007         ## $a2 = 00000007
/* 00818 80890028 00621824 */  and     $v1, $v1, $v0              
/* 0081C 8089002C 10600031 */  beq     $v1, $zero, .L808900F4     
/* 00820 80890030 3C0140C0 */  lui     $at, 0x40C0                ## $at = 40C00000
/* 00824 80890034 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 00828 80890038 C6080060 */  lwc1    $f8, 0x0060($s0)           ## 00000060
/* 0082C 8089003C C6040024 */  lwc1    $f4, 0x0024($s0)           ## 00000024
/* 00830 80890040 C606002C */  lwc1    $f6, 0x002C($s0)           ## 0000002C
/* 00834 80890044 4600403E */  c.le.s  $f8, $f0                   
/* 00838 80890048 E6000068 */  swc1    $f0, 0x0068($s0)           ## 00000068
/* 0083C 8089004C E6040168 */  swc1    $f4, 0x0168($s0)           ## 00000168
/* 00840 80890050 E6060170 */  swc1    $f6, 0x0170($s0)           ## 00000170
/* 00844 80890054 45020006 */  bc1fl   .L80890070                 
/* 00848 80890058 A600001C */  sh      $zero, 0x001C($s0)         ## 0000001C
/* 0084C 8089005C 8E0E0004 */  lw      $t6, 0x0004($s0)           ## 00000004
/* 00850 80890060 2401FFEF */  addiu   $at, $zero, 0xFFEF         ## $at = FFFFFFEF
/* 00854 80890064 01C17824 */  and     $t7, $t6, $at              
/* 00858 80890068 AE0F0004 */  sw      $t7, 0x0004($s0)           ## 00000004
/* 0085C 8089006C A600001C */  sh      $zero, 0x001C($s0)         ## 0000001C
.L80890070:
/* 00860 80890070 0C00B7D5 */  jal     func_8002DF54              
/* 00864 80890074 8FA4006C */  lw      $a0, 0x006C($sp)           
/* 00868 80890078 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 0086C 8089007C 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 00870 80890080 24052835 */  addiu   $a1, $zero, 0x2835         ## $a1 = 00002835
/* 00874 80890084 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 00878 80890088 44811000 */  mtc1    $at, $f2                   ## $f2 = 1.00
/* 0087C 8089008C 3C018089 */  lui     $at, %hi(D_80890538)       ## $at = 80890000
/* 00880 80890090 C4300538 */  lwc1    $f16, %lo(D_80890538)($at) 
/* 00884 80890094 C60A0024 */  lwc1    $f10, 0x0024($s0)          ## 00000024
/* 00888 80890098 3C198089 */  lui     $t9, %hi(func_8088FED0)    ## $t9 = 80890000
/* 0088C 8089009C 3C014382 */  lui     $at, 0x4382                ## $at = 43820000
/* 00890 808900A0 46105000 */  add.s   $f0, $f10, $f16            
/* 00894 808900A4 2739FED0 */  addiu   $t9, $t9, %lo(func_8088FED0) ## $t9 = 8088FED0
/* 00898 808900A8 46000005 */  abs.s   $f0, $f0                   
/* 0089C 808900AC 4602003C */  c.lt.s  $f0, $f2                   
/* 008A0 808900B0 00000000 */  nop
/* 008A4 808900B4 4500000D */  bc1f    .L808900EC                 
/* 008A8 808900B8 00000000 */  nop
/* 008AC 808900BC C612002C */  lwc1    $f18, 0x002C($s0)          ## 0000002C
/* 008B0 808900C0 44812000 */  mtc1    $at, $f4                   ## $f4 = 260.00
/* 008B4 808900C4 3C188089 */  lui     $t8, %hi(func_808903FC)    ## $t8 = 80890000
/* 008B8 808900C8 271803FC */  addiu   $t8, $t8, %lo(func_808903FC) ## $t8 = 808903FC
/* 008BC 808900CC 46049000 */  add.s   $f0, $f18, $f4             
/* 008C0 808900D0 46000005 */  abs.s   $f0, $f0                   
/* 008C4 808900D4 4602003C */  c.lt.s  $f0, $f2                   
/* 008C8 808900D8 00000000 */  nop
/* 008CC 808900DC 45000003 */  bc1f    .L808900EC                 
/* 008D0 808900E0 00000000 */  nop
/* 008D4 808900E4 10000096 */  beq     $zero, $zero, .L80890340   
/* 008D8 808900E8 AE180164 */  sw      $t8, 0x0164($s0)           ## 00000164
.L808900EC:
/* 008DC 808900EC 10000094 */  beq     $zero, $zero, .L80890340   
/* 008E0 808900F0 AE190164 */  sw      $t9, 0x0164($s0)           ## 00000164
.L808900F4:
/* 008E4 808900F4 44813000 */  mtc1    $at, $f6                   ## $f6 = 260.00
/* 008E8 808900F8 C6080068 */  lwc1    $f8, 0x0068($s0)           ## 00000068
/* 008EC 808900FC 4608303C */  c.lt.s  $f6, $f8                   
/* 008F0 80890100 00000000 */  nop
/* 008F4 80890104 4502008F */  bc1fl   .L80890344                 
/* 008F8 80890108 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 008FC 8089010C 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 00900 80890110 C60A0028 */  lwc1    $f10, 0x0028($s0)          ## 00000028
/* 00904 80890114 3C0142F0 */  lui     $at, 0x42F0                ## $at = 42F00000
/* 00908 80890118 460A003E */  c.le.s  $f0, $f10                  
/* 0090C 8089011C 00000000 */  nop
/* 00910 80890120 45020088 */  bc1fl   .L80890344                 
/* 00914 80890124 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00918 80890128 44816000 */  mtc1    $at, $f12                  ## $f12 = 120.00
/* 0091C 8089012C 0C00CFC8 */  jal     Math_Rand_CenteredFloat
              
/* 00920 80890130 00000000 */  nop
/* 00924 80890134 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 00928 80890138 E7A00048 */  swc1    $f0, 0x0048($sp)           
/* 0092C 8089013C 46000506 */  mov.s   $f20, $f0                  
/* 00930 80890140 0C01DE1C */  jal     Math_Sins
              ## sins?
/* 00934 80890144 86040158 */  lh      $a0, 0x0158($s0)           ## 00000158
/* 00938 80890148 3C013FC0 */  lui     $at, 0x3FC0                ## $at = 3FC00000
/* 0093C 8089014C 44818000 */  mtc1    $at, $f16                  ## $f16 = 1.50
/* 00940 80890150 00000000 */  nop
/* 00944 80890154 46148480 */  add.s   $f18, $f16, $f20           
/* 00948 80890158 46009107 */  neg.s   $f4, $f18                  
/* 0094C 8089015C 46040182 */  mul.s   $f6, $f0, $f4              
/* 00950 80890160 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 00954 80890164 E7A6004C */  swc1    $f6, 0x004C($sp)           
/* 00958 80890168 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 0095C 8089016C 44811000 */  mtc1    $at, $f2                   ## $f2 = 1.00
/* 00960 80890170 00000000 */  nop
/* 00964 80890174 46020200 */  add.s   $f8, $f0, $f2              
/* 00968 80890178 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 0096C 8089017C E7A80050 */  swc1    $f8, 0x0050($sp)           
/* 00970 80890180 46000506 */  mov.s   $f20, $f0                  
/* 00974 80890184 0C01DE0D */  jal     Math_Coss
              ## coss?
/* 00978 80890188 86040158 */  lh      $a0, 0x0158($s0)           ## 00000158
/* 0097C 8089018C 3C013FC0 */  lui     $at, 0x3FC0                ## $at = 3FC00000
/* 00980 80890190 44815000 */  mtc1    $at, $f10                  ## $f10 = 1.50
/* 00984 80890194 00000000 */  nop
/* 00988 80890198 46145400 */  add.s   $f16, $f10, $f20           
/* 0098C 8089019C 46008487 */  neg.s   $f18, $f16                 
/* 00990 808901A0 46120102 */  mul.s   $f4, $f0, $f18             
/* 00994 808901A4 E7A40054 */  swc1    $f4, 0x0054($sp)           
/* 00998 808901A8 0C01DE1C */  jal     Math_Sins
              ## sins?
/* 0099C 808901AC 86040158 */  lh      $a0, 0x0158($s0)           ## 00000158
/* 009A0 808901B0 46000506 */  mov.s   $f20, $f0                  
/* 009A4 808901B4 0C01DE0D */  jal     Math_Coss
              ## coss?
/* 009A8 808901B8 86040158 */  lh      $a0, 0x0158($s0)           ## 00000158
/* 009AC 808901BC 3C014270 */  lui     $at, 0x4270                ## $at = 42700000
/* 009B0 808901C0 44814000 */  mtc1    $at, $f8                   ## $f8 = 60.00
/* 009B4 808901C4 C7B20048 */  lwc1    $f18, 0x0048($sp)          
/* 009B8 808901C8 C6060024 */  lwc1    $f6, 0x0024($s0)           ## 00000024
/* 009BC 808901CC 46144282 */  mul.s   $f10, $f8, $f20            
/* 009C0 808901D0 460A3401 */  sub.s   $f16, $f6, $f10            
/* 009C4 808901D4 46120102 */  mul.s   $f4, $f0, $f18             
/* 009C8 808901D8 46048201 */  sub.s   $f8, $f16, $f4             
/* 009CC 808901DC E7A80058 */  swc1    $f8, 0x0058($sp)           
/* 009D0 808901E0 0C01DE0D */  jal     Math_Coss
              ## coss?
/* 009D4 808901E4 86040158 */  lh      $a0, 0x0158($s0)           ## 00000158
/* 009D8 808901E8 46000506 */  mov.s   $f20, $f0                  
/* 009DC 808901EC 0C01DE1C */  jal     Math_Sins
              ## sins?
/* 009E0 808901F0 86040158 */  lh      $a0, 0x0158($s0)           ## 00000158
/* 009E4 808901F4 3C014270 */  lui     $at, 0x4270                ## $at = 42700000
/* 009E8 808901F8 44815000 */  mtc1    $at, $f10                  ## $f10 = 60.00
/* 009EC 808901FC C7A40048 */  lwc1    $f4, 0x0048($sp)           
/* 009F0 80890200 C606002C */  lwc1    $f6, 0x002C($s0)           ## 0000002C
/* 009F4 80890204 46145482 */  mul.s   $f18, $f10, $f20           
/* 009F8 80890208 24040028 */  addiu   $a0, $zero, 0x0028         ## $a0 = 00000028
/* 009FC 8089020C 2405000F */  addiu   $a1, $zero, 0x000F         ## $a1 = 0000000F
/* 00A00 80890210 46040202 */  mul.s   $f8, $f0, $f4              
/* 00A04 80890214 46123401 */  sub.s   $f16, $f6, $f18            
/* 00A08 80890218 46104280 */  add.s   $f10, $f8, $f16            
/* 00A0C 8089021C E7AA0060 */  swc1    $f10, 0x0060($sp)          
/* 00A10 80890220 C6060028 */  lwc1    $f6, 0x0028($s0)           ## 00000028
/* 00A14 80890224 0C01DF64 */  jal     Math_Rand_S16Offset
              
/* 00A18 80890228 E7A6005C */  swc1    $f6, 0x005C($sp)           
/* 00A1C 8089022C 3C088089 */  lui     $t0, %hi(D_80890490)       ## $t0 = 80890000
/* 00A20 80890230 3C098089 */  lui     $t1, %hi(D_80890494)       ## $t1 = 80890000
/* 00A24 80890234 25290494 */  addiu   $t1, $t1, %lo(D_80890494)  ## $t1 = 80890494
/* 00A28 80890238 25080490 */  addiu   $t0, $t0, %lo(D_80890490)  ## $t0 = 80890490
/* 00A2C 8089023C 3C078089 */  lui     $a3, %hi(D_80890498)       ## $a3 = 80890000
/* 00A30 80890240 240A00FA */  addiu   $t2, $zero, 0x00FA         ## $t2 = 000000FA
/* 00A34 80890244 AFAA0018 */  sw      $t2, 0x0018($sp)           
/* 00A38 80890248 24E70498 */  addiu   $a3, $a3, %lo(D_80890498)  ## $a3 = 80890498
/* 00A3C 8089024C AFA80010 */  sw      $t0, 0x0010($sp)           
/* 00A40 80890250 AFA90014 */  sw      $t1, 0x0014($sp)           
/* 00A44 80890254 8FA4006C */  lw      $a0, 0x006C($sp)           
/* 00A48 80890258 27A50058 */  addiu   $a1, $sp, 0x0058           ## $a1 = FFFFFFF0
/* 00A4C 8089025C 27A6004C */  addiu   $a2, $sp, 0x004C           ## $a2 = FFFFFFE4
/* 00A50 80890260 0C00A0A7 */  jal     func_8002829C              
/* 00A54 80890264 AFA2001C */  sw      $v0, 0x001C($sp)           
/* 00A58 80890268 3C0142F0 */  lui     $at, 0x42F0                ## $at = 42F00000
/* 00A5C 8089026C 44816000 */  mtc1    $at, $f12                  ## $f12 = 120.00
/* 00A60 80890270 0C00CFC8 */  jal     Math_Rand_CenteredFloat
              
/* 00A64 80890274 00000000 */  nop
/* 00A68 80890278 E7A00048 */  swc1    $f0, 0x0048($sp)           
/* 00A6C 8089027C 0C01DE1C */  jal     Math_Sins
              ## sins?
/* 00A70 80890280 86040158 */  lh      $a0, 0x0158($s0)           ## 00000158
/* 00A74 80890284 46000506 */  mov.s   $f20, $f0                  
/* 00A78 80890288 0C01DE0D */  jal     Math_Coss
              ## coss?
/* 00A7C 8089028C 86040158 */  lh      $a0, 0x0158($s0)           ## 00000158
/* 00A80 80890290 3C014270 */  lui     $at, 0x4270                ## $at = 42700000
/* 00A84 80890294 44812000 */  mtc1    $at, $f4                   ## $f4 = 60.00
/* 00A88 80890298 C7AA0048 */  lwc1    $f10, 0x0048($sp)          
/* 00A8C 8089029C C6120024 */  lwc1    $f18, 0x0024($s0)          ## 00000024
/* 00A90 808902A0 46142202 */  mul.s   $f8, $f4, $f20             
/* 00A94 808902A4 46089401 */  sub.s   $f16, $f18, $f8            
/* 00A98 808902A8 460A0182 */  mul.s   $f6, $f0, $f10             
/* 00A9C 808902AC 46103100 */  add.s   $f4, $f6, $f16             
/* 00AA0 808902B0 E7A40058 */  swc1    $f4, 0x0058($sp)           
/* 00AA4 808902B4 0C01DE0D */  jal     Math_Coss
              ## coss?
/* 00AA8 808902B8 86040158 */  lh      $a0, 0x0158($s0)           ## 00000158
/* 00AAC 808902BC 46000506 */  mov.s   $f20, $f0                  
/* 00AB0 808902C0 0C01DE1C */  jal     Math_Sins
              ## sins?
/* 00AB4 808902C4 86040158 */  lh      $a0, 0x0158($s0)           ## 00000158
/* 00AB8 808902C8 3C014270 */  lui     $at, 0x4270                ## $at = 42700000
/* 00ABC 808902CC 44814000 */  mtc1    $at, $f8                   ## $f8 = 60.00
/* 00AC0 808902D0 C7B00048 */  lwc1    $f16, 0x0048($sp)          
/* 00AC4 808902D4 C612002C */  lwc1    $f18, 0x002C($s0)          ## 0000002C
/* 00AC8 808902D8 46144282 */  mul.s   $f10, $f8, $f20            
/* 00ACC 808902DC 24040028 */  addiu   $a0, $zero, 0x0028         ## $a0 = 00000028
/* 00AD0 808902E0 2405000F */  addiu   $a1, $zero, 0x000F         ## $a1 = 0000000F
/* 00AD4 808902E4 46100102 */  mul.s   $f4, $f0, $f16             
/* 00AD8 808902E8 460A9181 */  sub.s   $f6, $f18, $f10            
/* 00ADC 808902EC 46043201 */  sub.s   $f8, $f6, $f4              
/* 00AE0 808902F0 0C01DF64 */  jal     Math_Rand_S16Offset
              
/* 00AE4 808902F4 E7A80060 */  swc1    $f8, 0x0060($sp)           
/* 00AE8 808902F8 3C0B8089 */  lui     $t3, %hi(D_80890490)       ## $t3 = 80890000
/* 00AEC 808902FC 3C0C8089 */  lui     $t4, %hi(D_80890494)       ## $t4 = 80890000
/* 00AF0 80890300 258C0494 */  addiu   $t4, $t4, %lo(D_80890494)  ## $t4 = 80890494
/* 00AF4 80890304 256B0490 */  addiu   $t3, $t3, %lo(D_80890490)  ## $t3 = 80890490
/* 00AF8 80890308 3C078089 */  lui     $a3, %hi(D_80890498)       ## $a3 = 80890000
/* 00AFC 8089030C 240D00FA */  addiu   $t5, $zero, 0x00FA         ## $t5 = 000000FA
/* 00B00 80890310 AFAD0018 */  sw      $t5, 0x0018($sp)           
/* 00B04 80890314 24E70498 */  addiu   $a3, $a3, %lo(D_80890498)  ## $a3 = 80890498
/* 00B08 80890318 AFAB0010 */  sw      $t3, 0x0010($sp)           
/* 00B0C 8089031C AFAC0014 */  sw      $t4, 0x0014($sp)           
/* 00B10 80890320 8FA4006C */  lw      $a0, 0x006C($sp)           
/* 00B14 80890324 27A50058 */  addiu   $a1, $sp, 0x0058           ## $a1 = FFFFFFF0
/* 00B18 80890328 27A6004C */  addiu   $a2, $sp, 0x004C           ## $a2 = FFFFFFE4
/* 00B1C 8089032C 0C00A0A7 */  jal     func_8002829C              
/* 00B20 80890330 AFA2001C */  sw      $v0, 0x001C($sp)           
/* 00B24 80890334 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00B28 80890338 0C00BE5D */  jal     func_8002F974              
/* 00B2C 8089033C 240500DF */  addiu   $a1, $zero, 0x00DF         ## $a1 = 000000DF
.L80890340:
/* 00B30 80890340 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
.L80890344:
/* 00B34 80890344 0C223F52 */  jal     func_8088FD48              
/* 00B38 80890348 8FA5006C */  lw      $a1, 0x006C($sp)           
/* 00B3C 8089034C 8FBF0034 */  lw      $ra, 0x0034($sp)           
/* 00B40 80890350 D7B40028 */  ldc1    $f20, 0x0028($sp)          
/* 00B44 80890354 8FB00030 */  lw      $s0, 0x0030($sp)           
/* 00B48 80890358 03E00008 */  jr      $ra                        
/* 00B4C 8089035C 27BD0068 */  addiu   $sp, $sp, 0x0068           ## $sp = 00000000


