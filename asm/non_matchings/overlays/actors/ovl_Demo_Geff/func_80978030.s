.late_rodata

glabel jtbl_80978778
.word L80978064
.word func_809780E0
.word L8097816C
.word L80978064
.word func_809780E0
.word L8097816C
.word L80978064
.word func_809780E0
glabel D_80978798
 .word 0x3C23D70A
glabel D_8097879C
 .word 0x3E99999A

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


