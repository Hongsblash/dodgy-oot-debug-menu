glabel func_80A190C8
/* 01BB8 80A190C8 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 01BBC 80A190CC AFB00018 */  sw      $s0, 0x0018($sp)           
/* 01BC0 80A190D0 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 01BC4 80A190D4 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 01BC8 80A190D8 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 01BCC 80A190DC AFA50034 */  sw      $a1, 0x0034($sp)           
/* 01BD0 80A190E0 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 01BD4 80A190E4 AFA40020 */  sw      $a0, 0x0020($sp)           
/* 01BD8 80A190E8 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 01BDC 80A190EC 0C0295B2 */  jal     func_800A56C8              
/* 01BE0 80A190F0 24050000 */  addiu   $a1, $zero, 0x0000         ## $a1 = 00000000
/* 01BE4 80A190F4 14400004 */  bne     $v0, $zero, .L80A19108     
/* 01BE8 80A190F8 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 01BEC 80A190FC 0C0295B2 */  jal     func_800A56C8              
/* 01BF0 80A19100 3C054190 */  lui     $a1, 0x4190                ## $a1 = 41900000
/* 01BF4 80A19104 10400003 */  beq     $v0, $zero, .L80A19114     
.L80A19108:
/* 01BF8 80A19108 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01BFC 80A1910C 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 01C00 80A19110 24053931 */  addiu   $a1, $zero, 0x3931         ## $a1 = 00003931
.L80A19114:
/* 01C04 80A19114 96020088 */  lhu     $v0, 0x0088($s0)           ## 00000088
/* 01C08 80A19118 30420008 */  andi    $v0, $v0, 0x0008           ## $v0 = 00000000
/* 01C0C 80A1911C 50400008 */  beql    $v0, $zero, .L80A19140     
/* 01C10 80A19120 860F001C */  lh      $t7, 0x001C($s0)           ## 0000001C
/* 01C14 80A19124 860E007E */  lh      $t6, 0x007E($s0)           ## 0000007E
/* 01C18 80A19128 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01C1C 80A1912C 0C285E42 */  jal     func_80A17908              
/* 01C20 80A19130 A60E0196 */  sh      $t6, 0x0196($s0)           ## 00000196
/* 01C24 80A19134 1000003A */  beq     $zero, $zero, .L80A19220   
/* 01C28 80A19138 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 01C2C 80A1913C 860F001C */  lh      $t7, 0x001C($s0)           ## 0000001C
.L80A19140:
/* 01C30 80A19140 24010020 */  addiu   $at, $zero, 0x0020         ## $at = 00000020
/* 01C34 80A19144 260400B6 */  addiu   $a0, $s0, 0x00B6           ## $a0 = 000000B6
/* 01C38 80A19148 15E10028 */  bne     $t7, $at, .L80A191EC       
/* 01C3C 80A1914C 24060003 */  addiu   $a2, $zero, 0x0003         ## $a2 = 00000003
/* 01C40 80A19150 8E020118 */  lw      $v0, 0x0118($s0)           ## 00000118
/* 01C44 80A19154 24030040 */  addiu   $v1, $zero, 0x0040         ## $v1 = 00000040
/* 01C48 80A19158 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01C4C 80A1915C 8458001C */  lh      $t8, 0x001C($v0)           ## 0000001C
/* 01C50 80A19160 54780004 */  bnel    $v1, $t8, .L80A19174       
/* 01C54 80A19164 8E02011C */  lw      $v0, 0x011C($s0)           ## 0000011C
/* 01C58 80A19168 1000000A */  beq     $zero, $zero, .L80A19194   
/* 01C5C 80A1916C 00403825 */  or      $a3, $v0, $zero            ## $a3 = 00000000
/* 01C60 80A19170 8E02011C */  lw      $v0, 0x011C($s0)           ## 0000011C
.L80A19174:
/* 01C64 80A19174 24080010 */  addiu   $t0, $zero, 0x0010         ## $t0 = 00000010
/* 01C68 80A19178 8459001C */  lh      $t9, 0x001C($v0)           ## 0000001C
/* 01C6C 80A1917C 14790003 */  bne     $v1, $t9, .L80A1918C       
/* 01C70 80A19180 00000000 */  nop
/* 01C74 80A19184 10000003 */  beq     $zero, $zero, .L80A19194   
/* 01C78 80A19188 00403825 */  or      $a3, $v0, $zero            ## $a3 = 00000000
.L80A1918C:
/* 01C7C 80A1918C 10000023 */  beq     $zero, $zero, .L80A1921C   
/* 01C80 80A19190 A608001C */  sh      $t0, 0x001C($s0)           ## 0000001C
.L80A19194:
/* 01C84 80A19194 00E02825 */  or      $a1, $a3, $zero            ## $a1 = 00000000
/* 01C88 80A19198 0C00B69E */  jal     func_8002DA78              
/* 01C8C 80A1919C AFA7002C */  sw      $a3, 0x002C($sp)           
/* 01C90 80A191A0 00022C00 */  sll     $a1, $v0, 16               
/* 01C94 80A191A4 00052C03 */  sra     $a1, $a1, 16               
/* 01C98 80A191A8 260400B6 */  addiu   $a0, $s0, 0x00B6           ## $a0 = 000000B6
/* 01C9C 80A191AC 0C01DE2B */  jal     Math_ApproxUpdateScaledS
              
/* 01CA0 80A191B0 2406038E */  addiu   $a2, $zero, 0x038E         ## $a2 = 0000038E
/* 01CA4 80A191B4 8FA5002C */  lw      $a1, 0x002C($sp)           
/* 01CA8 80A191B8 0C00B6E3 */  jal     func_8002DB8C              
/* 01CAC 80A191BC 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01CB0 80A191C0 3C0142A0 */  lui     $at, 0x42A0                ## $at = 42A00000
/* 01CB4 80A191C4 44812000 */  mtc1    $at, $f4                   ## $f4 = 80.00
/* 01CB8 80A191C8 00000000 */  nop
/* 01CBC 80A191CC 4604003C */  c.lt.s  $f0, $f4                   
/* 01CC0 80A191D0 00000000 */  nop
/* 01CC4 80A191D4 45020012 */  bc1fl   .L80A19220                 
/* 01CC8 80A191D8 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 01CCC 80A191DC 0C285FB5 */  jal     func_80A17ED4              
/* 01CD0 80A191E0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01CD4 80A191E4 1000000E */  beq     $zero, $zero, .L80A19220   
/* 01CD8 80A191E8 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80A191EC:
/* 01CDC 80A191EC 8605008A */  lh      $a1, 0x008A($s0)           ## 0000008A
/* 01CE0 80A191F0 0C01E1EF */  jal     Math_SmoothScaleMaxS
              
/* 01CE4 80A191F4 2407071C */  addiu   $a3, $zero, 0x071C         ## $a3 = 0000071C
/* 01CE8 80A191F8 3C0142A0 */  lui     $at, 0x42A0                ## $at = 42A00000
/* 01CEC 80A191FC 44814000 */  mtc1    $at, $f8                   ## $f8 = 80.00
/* 01CF0 80A19200 C6060090 */  lwc1    $f6, 0x0090($s0)           ## 00000090
/* 01CF4 80A19204 4608303C */  c.lt.s  $f6, $f8                   
/* 01CF8 80A19208 00000000 */  nop
/* 01CFC 80A1920C 45020004 */  bc1fl   .L80A19220                 
/* 01D00 80A19210 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 01D04 80A19214 0C285FCF */  jal     func_80A17F3C              
/* 01D08 80A19218 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
.L80A1921C:
/* 01D0C 80A1921C 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80A19220:
/* 01D10 80A19220 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 01D14 80A19224 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000
/* 01D18 80A19228 03E00008 */  jr      $ra                        
/* 01D1C 80A1922C 00000000 */  nop


