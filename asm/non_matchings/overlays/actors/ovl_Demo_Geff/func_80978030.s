.late_rodata

glabel jtbl_80978778
.word L80978064
.word L809780E0
.word L8097816C
.word L80978064
.word L809780E0
.word L8097816C
.word L80978064
.word L809780E0
glabel D_80978798
 .float 0.01
glabel D_8097879C
 .float 0.3
glabel D_809787A0
 .float 0.01
glabel D_809787A4
 .float 0.15
glabel D_809787A8
 .float 0.29
glabel D_809787AC
 .float 0.12
glabel D_809787B0
 .float 0.01
glabel D_809787B4
 .float 0.1
glabel D_809787B8
 .float 0.15
glabel D_809787BC
 .float 0.2

.text

glabel func_80978030
/* 001F0 80978030 AFA50004 */  sw      $a1, 0x0004($sp)           
/* 001F4 80978034 948E001C */  lhu     $t6, 0x001C($a0)           ## 0000001C
/* 001F8 80978038 24020001 */  addiu   $v0, $zero, 0x0001         ## $v0 = 00000001
/* 001FC 8097803C AC82014C */  sw      $v0, 0x014C($a0)           ## 0000014C
/* 00200 80978040 2DC10008 */  sltiu   $at, $t6, 0x0008           
/* 00204 80978044 10200049 */  beq     $at, $zero, .L8097816C     
/* 00208 80978048 AC820150 */  sw      $v0, 0x0150($a0)           ## 00000150
/* 0020C 8097804C 000E7080 */  sll     $t6, $t6,  2               
/* 00210 80978050 3C018098 */  lui     $at, %hi(jtbl_80978778)       ## $at = 80980000
/* 00214 80978054 002E0821 */  addu    $at, $at, $t6              
/* 00218 80978058 8C2E8778 */  lw      $t6, %lo(jtbl_80978778)($at)  
/* 0021C 8097805C 01C00008 */  jr      $t6                        
/* 00220 80978060 00000000 */  nop
glabel L80978064
.L80978064:
/* 00224 80978064 3C038016 */  lui     $v1, 0x8016                ## $v1 = 80160000
/* 00228 80978068 2463FA90 */  addiu   $v1, $v1, 0xFA90           ## $v1 = 8015FA90
/* 0022C 8097806C 8C6F0000 */  lw      $t7, 0x0000($v1)           ## 8015FA90
/* 00230 80978070 3C018098 */  lui     $at, %hi(D_80978798)       ## $at = 80980000
/* 00234 80978074 C4208798 */  lwc1    $f0, %lo(D_80978798)($at)  
/* 00238 80978078 85F81462 */  lh      $t8, 0x1462($t7)           ## 00001462
/* 0023C 8097807C 3C018098 */  lui     $at, %hi(D_8097879C)       ## $at = 80980000
/* 00240 80978080 C422879C */  lwc1    $f2, %lo(D_8097879C)($at)  
/* 00244 80978084 44982000 */  mtc1    $t8, $f4                   ## $f4 = 0.00
/* 00248 80978088 24820050 */  addiu   $v0, $a0, 0x0050           ## $v0 = 00000050
/* 0024C 8097808C 468021A0 */  cvt.s.w $f6, $f4                   
/* 00250 80978090 46003202 */  mul.s   $f8, $f6, $f0              
/* 00254 80978094 46024280 */  add.s   $f10, $f8, $f2             
/* 00258 80978098 E44A0000 */  swc1    $f10, 0x0000($v0)          ## 00000050
/* 0025C 8097809C 8C790000 */  lw      $t9, 0x0000($v1)           ## 8015FA90
/* 00260 809780A0 87281464 */  lh      $t0, 0x1464($t9)           ## 00001464
/* 00264 809780A4 44888000 */  mtc1    $t0, $f16                  ## $f16 = 0.00
/* 00268 809780A8 00000000 */  nop
/* 0026C 809780AC 468084A0 */  cvt.s.w $f18, $f16                 
/* 00270 809780B0 46009102 */  mul.s   $f4, $f18, $f0             
/* 00274 809780B4 46022180 */  add.s   $f6, $f4, $f2              
/* 00278 809780B8 E4460004 */  swc1    $f6, 0x0004($v0)           ## 00000054
/* 0027C 809780BC 8C690000 */  lw      $t1, 0x0000($v1)           ## 8015FA90
/* 00280 809780C0 852A1466 */  lh      $t2, 0x1466($t1)           ## 00001466
/* 00284 809780C4 448A4000 */  mtc1    $t2, $f8                   ## $f8 = 0.00
/* 00288 809780C8 00000000 */  nop
/* 0028C 809780CC 468042A0 */  cvt.s.w $f10, $f8                  
/* 00290 809780D0 46005402 */  mul.s   $f16, $f10, $f0            
/* 00294 809780D4 46028480 */  add.s   $f18, $f16, $f2            
/* 00298 809780D8 03E00008 */  jr      $ra                        
/* 0029C 809780DC E4520008 */  swc1    $f18, 0x0008($v0)          ## 00000058
glabel L809780E0
.L809780E0:
/* 002A0 809780E0 3C038016 */  lui     $v1, 0x8016                ## $v1 = 80160000
/* 002A4 809780E4 2463FA90 */  addiu   $v1, $v1, 0xFA90           ## $v1 = 8015FA90
/* 002A8 809780E8 8C6B0000 */  lw      $t3, 0x0000($v1)           ## 8015FA90
/* 002AC 809780EC 3C018098 */  lui     $at, %hi(D_809787A0)       ## $at = 80980000
/* 002B0 809780F0 C42087A0 */  lwc1    $f0, %lo(D_809787A0)($at)  
/* 002B4 809780F4 856C1468 */  lh      $t4, 0x1468($t3)           ## 00001468
/* 002B8 809780F8 3C018098 */  lui     $at, %hi(D_809787A4)       ## $at = 80980000
/* 002BC 809780FC C42A87A4 */  lwc1    $f10, %lo(D_809787A4)($at) 
/* 002C0 80978100 448C2000 */  mtc1    $t4, $f4                   ## $f4 = 0.00
/* 002C4 80978104 24820050 */  addiu   $v0, $a0, 0x0050           ## $v0 = 00000050
/* 002C8 80978108 3C018098 */  lui     $at, %hi(D_809787A8)       ## $at = 80980000
/* 002CC 8097810C 468021A0 */  cvt.s.w $f6, $f4                   
/* 002D0 80978110 46003202 */  mul.s   $f8, $f6, $f0              
/* 002D4 80978114 460A4400 */  add.s   $f16, $f8, $f10            
/* 002D8 80978118 E4500000 */  swc1    $f16, 0x0000($v0)          ## 00000050
/* 002DC 8097811C 8C6D0000 */  lw      $t5, 0x0000($v1)           ## 8015FA90
/* 002E0 80978120 C42887A8 */  lwc1    $f8, %lo(D_809787A8)($at)  
/* 002E4 80978124 3C018098 */  lui     $at, %hi(D_809787AC)       ## $at = 80980000
/* 002E8 80978128 85AE146A */  lh      $t6, 0x146A($t5)           ## 0000146A
/* 002EC 8097812C 448E9000 */  mtc1    $t6, $f18                  ## $f18 = 0.00
/* 002F0 80978130 00000000 */  nop
/* 002F4 80978134 46809120 */  cvt.s.w $f4, $f18                  
/* 002F8 80978138 46002182 */  mul.s   $f6, $f4, $f0              
/* 002FC 8097813C 46083280 */  add.s   $f10, $f6, $f8             
/* 00300 80978140 E44A0004 */  swc1    $f10, 0x0004($v0)          ## 00000054
/* 00304 80978144 8C6F0000 */  lw      $t7, 0x0000($v1)           ## 8015FA90
/* 00308 80978148 C42687AC */  lwc1    $f6, %lo(D_809787AC)($at)  
/* 0030C 8097814C 85F8146C */  lh      $t8, 0x146C($t7)           ## 0000146C
/* 00310 80978150 44988000 */  mtc1    $t8, $f16                  ## $f16 = 0.00
/* 00314 80978154 00000000 */  nop
/* 00318 80978158 468084A0 */  cvt.s.w $f18, $f16                 
/* 0031C 8097815C 46009102 */  mul.s   $f4, $f18, $f0             
/* 00320 80978160 46062200 */  add.s   $f8, $f4, $f6              
/* 00324 80978164 03E00008 */  jr      $ra                        
/* 00328 80978168 E4480008 */  swc1    $f8, 0x0008($v0)           ## 00000058
glabel L8097816C
.L8097816C:
/* 0032C 8097816C 3C038016 */  lui     $v1, 0x8016                ## $v1 = 80160000
/* 00330 80978170 2463FA90 */  addiu   $v1, $v1, 0xFA90           ## $v1 = 8015FA90
/* 00334 80978174 8C790000 */  lw      $t9, 0x0000($v1)           ## 8015FA90
/* 00338 80978178 3C018098 */  lui     $at, %hi(D_809787B0)       ## $at = 80980000
/* 0033C 8097817C C42087B0 */  lwc1    $f0, %lo(D_809787B0)($at)  
/* 00340 80978180 8728146E */  lh      $t0, 0x146E($t9)           ## 0000146E
/* 00344 80978184 3C018098 */  lui     $at, %hi(D_809787B4)       ## $at = 80980000
/* 00348 80978188 C42487B4 */  lwc1    $f4, %lo(D_809787B4)($at)  
/* 0034C 8097818C 44885000 */  mtc1    $t0, $f10                  ## $f10 = 0.00
/* 00350 80978190 24820050 */  addiu   $v0, $a0, 0x0050           ## $v0 = 00000050
/* 00354 80978194 3C018098 */  lui     $at, %hi(D_809787B8)       ## $at = 80980000
/* 00358 80978198 46805420 */  cvt.s.w $f16, $f10                 
/* 0035C 8097819C 46008482 */  mul.s   $f18, $f16, $f0            
/* 00360 809781A0 46049180 */  add.s   $f6, $f18, $f4             
/* 00364 809781A4 E4460000 */  swc1    $f6, 0x0000($v0)           ## 00000050
/* 00368 809781A8 8C690000 */  lw      $t1, 0x0000($v1)           ## 8015FA90
/* 0036C 809781AC C43287B8 */  lwc1    $f18, %lo(D_809787B8)($at) 
/* 00370 809781B0 3C018098 */  lui     $at, %hi(D_809787BC)       ## $at = 80980000
/* 00374 809781B4 852A1470 */  lh      $t2, 0x1470($t1)           ## 00001470
/* 00378 809781B8 448A4000 */  mtc1    $t2, $f8                   ## $f8 = 0.00
/* 0037C 809781BC 00000000 */  nop
/* 00380 809781C0 468042A0 */  cvt.s.w $f10, $f8                  
/* 00384 809781C4 46005402 */  mul.s   $f16, $f10, $f0            
/* 00388 809781C8 46128100 */  add.s   $f4, $f16, $f18            
/* 0038C 809781CC E4440004 */  swc1    $f4, 0x0004($v0)           ## 00000054
/* 00390 809781D0 8C6B0000 */  lw      $t3, 0x0000($v1)           ## 8015FA90
/* 00394 809781D4 C43087BC */  lwc1    $f16, %lo(D_809787BC)($at) 
/* 00398 809781D8 856C1472 */  lh      $t4, 0x1472($t3)           ## 00001472
/* 0039C 809781DC 448C3000 */  mtc1    $t4, $f6                   ## $f6 = 0.00
/* 003A0 809781E0 00000000 */  nop
/* 003A4 809781E4 46803220 */  cvt.s.w $f8, $f6                   
/* 003A8 809781E8 46004282 */  mul.s   $f10, $f8, $f0             
/* 003AC 809781EC 46105480 */  add.s   $f18, $f10, $f16           
/* 003B0 809781F0 E4520008 */  swc1    $f18, 0x0008($v0)          ## 00000058
/* 003B4 809781F4 03E00008 */  jr      $ra                        
/* 003B8 809781F8 00000000 */  nop
