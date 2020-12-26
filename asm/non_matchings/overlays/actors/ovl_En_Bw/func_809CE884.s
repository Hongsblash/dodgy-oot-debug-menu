.late_rodata
glabel D_809D15A4
    .float 0.001

glabel D_809D15A8
 .word 0x3BE56041
glabel D_809D15AC
 .word 0x3C54FDF4
glabel D_809D15B0
    .float 0.001

glabel D_809D15B4
 .word 0x3C54FDF4
glabel D_809D15B8
 .word 0x3BE56041
glabel D_809D15BC
    .float 0.001

glabel D_809D15C0
 .word 0x3BE56041
glabel D_809D15C4
 .word 0x3C54FDF4

.text
glabel func_809CE884
/* 001A4 809CE884 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 001A8 809CE888 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 001AC 809CE88C AFB00020 */  sw      $s0, 0x0020($sp)           
/* 001B0 809CE890 AFA5002C */  sw      $a1, 0x002C($sp)           
/* 001B4 809CE894 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 001B8 809CE898 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 001BC 809CE89C 44050000 */  mfc1    $a1, $f0                   
/* 001C0 809CE8A0 24840068 */  addiu   $a0, $a0, 0x0068           ## $a0 = 00000068
/* 001C4 809CE8A4 3C063F80 */  lui     $a2, 0x3F80                ## $a2 = 3F800000
/* 001C8 809CE8A8 3C073F00 */  lui     $a3, 0x3F00                ## $a3 = 3F000000
/* 001CC 809CE8AC 0C01E0C4 */  jal     Math_SmoothStepToF
              
/* 001D0 809CE8B0 E7A00010 */  swc1    $f0, 0x0010($sp)           
/* 001D4 809CE8B4 860E0222 */  lh      $t6, 0x0222($s0)           ## 00000222
/* 001D8 809CE8B8 3C01809D */  lui     $at, %hi(D_809D15A4)       ## $at = 809D0000
/* 001DC 809CE8BC 25CFFF06 */  addiu   $t7, $t6, 0xFF06           ## $t7 = FFFFFF06
/* 001E0 809CE8C0 A60F0222 */  sh      $t7, 0x0222($s0)           ## 00000222
/* 001E4 809CE8C4 86180222 */  lh      $t8, 0x0222($s0)           ## 00000222
/* 001E8 809CE8C8 C42815A4 */  lwc1    $f8, %lo(D_809D15A4)($at)  
/* 001EC 809CE8CC 44982000 */  mtc1    $t8, $f4                   ## $f4 = 0.00
/* 001F0 809CE8D0 00000000 */  nop
/* 001F4 809CE8D4 468021A0 */  cvt.s.w $f6, $f4                   
/* 001F8 809CE8D8 46083302 */  mul.s   $f12, $f6, $f8             
/* 001FC 809CE8DC 0C0329C8 */  jal     Math_SinF              
/* 00200 809CE8E0 00000000 */  nop
/* 00204 809CE8E4 3C01809D */  lui     $at, %hi(D_809D15A8)       ## $at = 809D0000
/* 00208 809CE8E8 C42A15A8 */  lwc1    $f10, %lo(D_809D15A8)($at) 
/* 0020C 809CE8EC 86190222 */  lh      $t9, 0x0222($s0)           ## 00000222
/* 00210 809CE8F0 3C01809D */  lui     $at, %hi(D_809D15AC)       ## $at = 809D0000
/* 00214 809CE8F4 460A0402 */  mul.s   $f16, $f0, $f10            
/* 00218 809CE8F8 C43215AC */  lwc1    $f18, %lo(D_809D15AC)($at) 
/* 0021C 809CE8FC 44993000 */  mtc1    $t9, $f6                   ## $f6 = 0.00
/* 00220 809CE900 3C01809D */  lui     $at, %hi(D_809D15B0)       ## $at = 809D0000
/* 00224 809CE904 46803220 */  cvt.s.w $f8, $f6                   
/* 00228 809CE908 46128100 */  add.s   $f4, $f16, $f18            
/* 0022C 809CE90C E6040050 */  swc1    $f4, 0x0050($s0)           ## 00000050
/* 00230 809CE910 C42A15B0 */  lwc1    $f10, %lo(D_809D15B0)($at) 
/* 00234 809CE914 460A4302 */  mul.s   $f12, $f8, $f10            
/* 00238 809CE918 0C0329C8 */  jal     Math_SinF              
/* 0023C 809CE91C 00000000 */  nop
/* 00240 809CE920 3C01809D */  lui     $at, %hi(D_809D15B4)       ## $at = 809D0000
/* 00244 809CE924 C43015B4 */  lwc1    $f16, %lo(D_809D15B4)($at) 
/* 00248 809CE928 3C01809D */  lui     $at, %hi(D_809D15B8)       ## $at = 809D0000
/* 0024C 809CE92C C43215B8 */  lwc1    $f18, %lo(D_809D15B8)($at) 
/* 00250 809CE930 86080222 */  lh      $t0, 0x0222($s0)           ## 00000222
/* 00254 809CE934 3C01809D */  lui     $at, %hi(D_809D15BC)       ## $at = 809D0000
/* 00258 809CE938 46120102 */  mul.s   $f4, $f0, $f18             
/* 0025C 809CE93C 44884000 */  mtc1    $t0, $f8                   ## $f8 = 0.00
/* 00260 809CE940 00000000 */  nop
/* 00264 809CE944 468042A0 */  cvt.s.w $f10, $f8                  
/* 00268 809CE948 46048181 */  sub.s   $f6, $f16, $f4             
/* 0026C 809CE94C E6060054 */  swc1    $f6, 0x0054($s0)           ## 00000054
/* 00270 809CE950 C43215BC */  lwc1    $f18, %lo(D_809D15BC)($at) 
/* 00274 809CE954 46125302 */  mul.s   $f12, $f10, $f18           
/* 00278 809CE958 0C0329C8 */  jal     Math_SinF              
/* 0027C 809CE95C 00000000 */  nop
/* 00280 809CE960 3C01809D */  lui     $at, %hi(D_809D15C0)       ## $at = 809D0000
/* 00284 809CE964 C43015C0 */  lwc1    $f16, %lo(D_809D15C0)($at) 
/* 00288 809CE968 3C01809D */  lui     $at, %hi(D_809D15C4)       ## $at = 809D0000
/* 0028C 809CE96C C42615C4 */  lwc1    $f6, %lo(D_809D15C4)($at)  
/* 00290 809CE970 46100102 */  mul.s   $f4, $f0, $f16             
/* 00294 809CE974 86090222 */  lh      $t1, 0x0222($s0)           ## 00000222
/* 00298 809CE978 46062200 */  add.s   $f8, $f4, $f6              
/* 0029C 809CE97C 15200005 */  bne     $t1, $zero, .L809CE994     
/* 002A0 809CE980 E6080058 */  swc1    $f8, 0x0058($s0)           ## 00000058
/* 002A4 809CE984 860A00B6 */  lh      $t2, 0x00B6($s0)           ## 000000B6
/* 002A8 809CE988 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 002AC 809CE98C 0C273A6A */  jal     func_809CE9A8              
/* 002B0 809CE990 A60A0032 */  sh      $t2, 0x0032($s0)           ## 00000032
.L809CE994:
/* 002B4 809CE994 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 002B8 809CE998 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 002BC 809CE99C 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 002C0 809CE9A0 03E00008 */  jr      $ra                        
/* 002C4 809CE9A4 00000000 */  nop
