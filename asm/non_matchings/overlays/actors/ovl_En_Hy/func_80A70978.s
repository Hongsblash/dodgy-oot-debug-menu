glabel func_80A70978
/* 013C8 80A70978 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 013CC 80A7097C AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 013D0 80A70980 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 013D4 80A70984 AFA5003C */  sw      $a1, 0x003C($sp)           
/* 013D8 80A70988 848F001C */  lh      $t7, 0x001C($a0)           ## 0000001C
/* 013DC 80A7098C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 013E0 80A70990 8CA21C44 */  lw      $v0, 0x1C44($a1)           ## 00001C44
/* 013E4 80A70994 31F8007F */  andi    $t8, $t7, 0x007F           ## $t8 = 00000000
/* 013E8 80A70998 2F010013 */  sltiu   $at, $t8, 0x0013           
/* 013EC 80A7099C 10200018 */  beq     $at, $zero, .L80A70A00     
/* 013F0 80A709A0 0018C080 */  sll     $t8, $t8,  2               
/* 013F4 80A709A4 3C0180A7 */  lui     $at, %hi(jtbl_80A72980)       ## $at = 80A70000
/* 013F8 80A709A8 00380821 */  addu    $at, $at, $t8              
/* 013FC 80A709AC 8C382980 */  lw      $t8, %lo(jtbl_80A72980)($at)  
/* 01400 80A709B0 03000008 */  jr      $t8                        
/* 01404 80A709B4 00000000 */  nop
glabel L80A709B8
/* 01408 80A709B8 861901E8 */  lh      $t9, 0x01E8($s0)           ## 000001E8
/* 0140C 80A709BC 17200003 */  bne     $t9, $zero, .L80A709CC     
/* 01410 80A709C0 00000000 */  nop
/* 01414 80A709C4 1000000F */  beq     $zero, $zero, .L80A70A04   
/* 01418 80A709C8 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
.L80A709CC:
/* 0141C 80A709CC 1000000D */  beq     $zero, $zero, .L80A70A04   
/* 01420 80A709D0 24070002 */  addiu   $a3, $zero, 0x0002         ## $a3 = 00000002
glabel L80A709D4
/* 01424 80A709D4 1000000B */  beq     $zero, $zero, .L80A70A04   
/* 01428 80A709D8 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
glabel L80A709DC
/* 0142C 80A709DC 10000009 */  beq     $zero, $zero, .L80A70A04   
/* 01430 80A709E0 24070004 */  addiu   $a3, $zero, 0x0004         ## $a3 = 00000004
glabel L80A709E4
/* 01434 80A709E4 860801E8 */  lh      $t0, 0x01E8($s0)           ## 000001E8
/* 01438 80A709E8 15000003 */  bne     $t0, $zero, .L80A709F8     
/* 0143C 80A709EC 00000000 */  nop
/* 01440 80A709F0 10000004 */  beq     $zero, $zero, .L80A70A04   
/* 01444 80A709F4 24070002 */  addiu   $a3, $zero, 0x0002         ## $a3 = 00000002
.L80A709F8:
/* 01448 80A709F8 10000002 */  beq     $zero, $zero, .L80A70A04   
/* 0144C 80A709FC 24070004 */  addiu   $a3, $zero, 0x0004         ## $a3 = 00000004
glabel L80A70A00
.L80A70A00:
/* 01450 80A70A00 24070002 */  addiu   $a3, $zero, 0x0002         ## $a3 = 00000002
.L80A70A04:
/* 01454 80A70A04 8C4A0024 */  lw      $t2, 0x0024($v0)           ## 00000024
/* 01458 80A70A08 3C0B8016 */  lui     $t3, 0x8016                ## $t3 = 80160000
/* 0145C 80A70A0C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01460 80A70A10 AE0A0200 */  sw      $t2, 0x0200($s0)           ## 00000200
/* 01464 80A70A14 8C490028 */  lw      $t1, 0x0028($v0)           ## 00000028
/* 01468 80A70A18 260501E8 */  addiu   $a1, $s0, 0x01E8           ## $a1 = 000001E8
/* 0146C 80A70A1C AE090204 */  sw      $t1, 0x0204($s0)           ## 00000204
/* 01470 80A70A20 8C4A002C */  lw      $t2, 0x002C($v0)           ## 0000002C
/* 01474 80A70A24 AE0A0208 */  sw      $t2, 0x0208($s0)           ## 00000208
/* 01478 80A70A28 8D6BE664 */  lw      $t3, -0x199C($t3)          ## 8015E664
/* 0147C 80A70A2C 5560000D */  bnel    $t3, $zero, .L80A70A64     
/* 01480 80A70A30 8618001C */  lh      $t8, 0x001C($s0)           ## 0000001C
/* 01484 80A70A34 860C001C */  lh      $t4, 0x001C($s0)           ## 0000001C
/* 01488 80A70A38 3C0F80A7 */  lui     $t7, %hi(D_80A724A8)       ## $t7 = 80A70000
/* 0148C 80A70A3C 25EF24A8 */  addiu   $t7, $t7, %lo(D_80A724A8)  ## $t7 = 80A724A8
/* 01490 80A70A40 318D007F */  andi    $t5, $t4, 0x007F           ## $t5 = 00000000
/* 01494 80A70A44 000D7080 */  sll     $t6, $t5,  2               
/* 01498 80A70A48 01CD7023 */  subu    $t6, $t6, $t5              
/* 0149C 80A70A4C 000E7080 */  sll     $t6, $t6,  2               
/* 014A0 80A70A50 01CF1021 */  addu    $v0, $t6, $t7              
/* 014A4 80A70A54 C4440008 */  lwc1    $f4, 0x0008($v0)           ## 00000008
/* 014A8 80A70A58 1000000B */  beq     $zero, $zero, .L80A70A88   
/* 014AC 80A70A5C E60401FC */  swc1    $f4, 0x01FC($s0)           ## 000001FC
/* 014B0 80A70A60 8618001C */  lh      $t8, 0x001C($s0)           ## 0000001C
.L80A70A64:
/* 014B4 80A70A64 3C0980A7 */  lui     $t1, %hi(D_80A724A8)       ## $t1 = 80A70000
/* 014B8 80A70A68 252924A8 */  addiu   $t1, $t1, %lo(D_80A724A8)  ## $t1 = 80A724A8
/* 014BC 80A70A6C 3319007F */  andi    $t9, $t8, 0x007F           ## $t9 = 00000000
/* 014C0 80A70A70 00194080 */  sll     $t0, $t9,  2               
/* 014C4 80A70A74 01194023 */  subu    $t0, $t0, $t9              
/* 014C8 80A70A78 00084080 */  sll     $t0, $t0,  2               
/* 014CC 80A70A7C 01091021 */  addu    $v0, $t0, $t1              
/* 014D0 80A70A80 C4460004 */  lwc1    $f6, 0x0004($v0)           ## 00000004
/* 014D4 80A70A84 E60601FC */  swc1    $f6, 0x01FC($s0)           ## 000001FC
.L80A70A88:
/* 014D8 80A70A88 90460000 */  lbu     $a2, 0x0000($v0)           ## 00000000
/* 014DC 80A70A8C 0C00D285 */  jal     func_80034A14              
/* 014E0 80A70A90 AFA50028 */  sw      $a1, 0x0028($sp)           
/* 014E4 80A70A94 3C0A80A7 */  lui     $t2, %hi(func_80A6F810)    ## $t2 = 80A70000
/* 014E8 80A70A98 3C0B80A7 */  lui     $t3, %hi(func_80A70058)    ## $t3 = 80A70000
/* 014EC 80A70A9C 256B0058 */  addiu   $t3, $t3, %lo(func_80A70058) ## $t3 = 80A70058
/* 014F0 80A70AA0 254AF810 */  addiu   $t2, $t2, %lo(func_80A6F810) ## $t2 = 80A6F810
/* 014F4 80A70AA4 8E07025C */  lw      $a3, 0x025C($s0)           ## 0000025C
/* 014F8 80A70AA8 AFAB0014 */  sw      $t3, 0x0014($sp)           
/* 014FC 80A70AAC AFAA0010 */  sw      $t2, 0x0010($sp)           
/* 01500 80A70AB0 8FA4003C */  lw      $a0, 0x003C($sp)           
/* 01504 80A70AB4 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 01508 80A70AB8 0C00D0F3 */  jal     func_800343CC              
/* 0150C 80A70ABC 8FA60028 */  lw      $a2, 0x0028($sp)           
/* 01510 80A70AC0 10400003 */  beq     $v0, $zero, .L80A70AD0     
/* 01514 80A70AC4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01518 80A70AC8 0C29C20D */  jal     func_80A70834              
/* 0151C 80A70ACC 8FA5003C */  lw      $a1, 0x003C($sp)           
.L80A70AD0:
/* 01520 80A70AD0 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 01524 80A70AD4 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 01528 80A70AD8 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 0152C 80A70ADC 03E00008 */  jr      $ra                        
/* 01530 80A70AE0 00000000 */  nop


