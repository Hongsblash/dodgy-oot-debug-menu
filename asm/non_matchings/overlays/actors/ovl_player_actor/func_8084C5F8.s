glabel func_8084C5F8
/* 1A3E8 8084C5F8 27BDFFC0 */  addiu   $sp, $sp, 0xFFC0           ## $sp = FFFFFFC0
/* 1A3EC 8084C5FC AFBF001C */  sw      $ra, 0x001C($sp)           
/* 1A3F0 8084C600 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 1A3F4 8084C604 AFA50044 */  sw      $a1, 0x0044($sp)           
/* 1A3F8 8084C608 8C8E0680 */  lw      $t6, 0x0680($a0)           ## 00000680
/* 1A3FC 8084C60C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 1A400 8084C610 260601B4 */  addiu   $a2, $s0, 0x01B4           ## $a2 = 000001B4
/* 1A404 8084C614 35CF0040 */  ori     $t7, $t6, 0x0040           ## $t7 = 00000040
/* 1A408 8084C618 AC8F0680 */  sw      $t7, 0x0680($a0)           ## 00000680
/* 1A40C 8084C61C AFA60020 */  sw      $a2, 0x0020($sp)           
/* 1A410 8084C620 8FA40044 */  lw      $a0, 0x0044($sp)           
/* 1A414 8084C624 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 1A418 8084C628 0C20DD28 */  jal     func_808374A0              
/* 1A41C 8084C62C 3C074080 */  lui     $a3, 0x4080                ## $a3 = 40800000
/* 1A420 8084C630 14400007 */  bne     $v0, $zero, .L8084C650     
/* 1A424 8084C634 00000000 */  nop
/* 1A428 8084C638 8E18067C */  lw      $t8, 0x067C($s0)           ## 0000067C
/* 1A42C 8084C63C 3C01FFDF */  lui     $at, 0xFFDF                ## $at = FFDF0000
/* 1A430 8084C640 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = FFDFFFFF
/* 1A434 8084C644 0301C824 */  and     $t9, $t8, $at              
/* 1A438 8084C648 10000040 */  beq     $zero, $zero, .L8084C74C   
/* 1A43C 8084C64C AE19067C */  sw      $t9, 0x067C($s0)           ## 0000067C
.L8084C650:
/* 1A440 8084C650 1C400004 */  bgtz    $v0, .L8084C664            
/* 1A444 8084C654 8FA40044 */  lw      $a0, 0x0044($sp)           
/* 1A448 8084C658 0C028EF0 */  jal     func_800A3BC0              
/* 1A44C 8084C65C 8FA50020 */  lw      $a1, 0x0020($sp)           
/* 1A450 8084C660 10400009 */  beq     $v0, $zero, .L8084C688     
.L8084C664:
/* 1A454 8084C664 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 1A458 8084C668 0C20F03A */  jal     func_8083C0E8              
/* 1A45C 8084C66C 8FA50044 */  lw      $a1, 0x0044($sp)           
/* 1A460 8084C670 8E08067C */  lw      $t0, 0x067C($s0)           ## 0000067C
/* 1A464 8084C674 3C01FFDF */  lui     $at, 0xFFDF                ## $at = FFDF0000
/* 1A468 8084C678 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = FFDFFFFF
/* 1A46C 8084C67C 01014824 */  and     $t1, $t0, $at              
/* 1A470 8084C680 10000032 */  beq     $zero, $zero, .L8084C74C   
/* 1A474 8084C684 AE09067C */  sw      $t1, 0x067C($s0)           ## 0000067C
.L8084C688:
/* 1A478 8084C688 860A0850 */  lh      $t2, 0x0850($s0)           ## 00000850
/* 1A47C 8084C68C 3C038085 */  lui     $v1, %hi(D_80854898)       ## $v1 = 80850000
/* 1A480 8084C690 24634898 */  addiu   $v1, $v1, %lo(D_80854898)  ## $v1 = 80854898
/* 1A484 8084C694 11400006 */  beq     $t2, $zero, .L8084C6B0     
/* 1A488 8084C698 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 1A48C 8084C69C 3C058085 */  lui     $a1, %hi(D_808548A8)       ## $a1 = 80850000
/* 1A490 8084C6A0 0C20CA49 */  jal     func_80832924              
/* 1A494 8084C6A4 24A548A8 */  addiu   $a1, $a1, %lo(D_808548A8)  ## $a1 = 808548A8
/* 1A498 8084C6A8 3C038085 */  lui     $v1, %hi(D_808548A0)       ## $v1 = 80850000
/* 1A49C 8084C6AC 246348A0 */  addiu   $v1, $v1, %lo(D_808548A0)  ## $v1 = 808548A0
.L8084C6B0:
/* 1A4A0 8084C6B0 8C650000 */  lw      $a1, 0x0000($v1)           ## 808548A0
/* 1A4A4 8084C6B4 AFA30038 */  sw      $v1, 0x0038($sp)           
/* 1A4A8 8084C6B8 0C02914C */  jal     func_800A4530              
/* 1A4AC 8084C6BC 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 1A4B0 8084C6C0 14400006 */  bne     $v0, $zero, .L8084C6DC     
/* 1A4B4 8084C6C4 8FA30038 */  lw      $v1, 0x0038($sp)           
/* 1A4B8 8084C6C8 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 1A4BC 8084C6CC 0C02914C */  jal     func_800A4530              
/* 1A4C0 8084C6D0 8C650004 */  lw      $a1, 0x0004($v1)           ## 00000004
/* 1A4C4 8084C6D4 5040001E */  beql    $v0, $zero, .L8084C750     
/* 1A4C8 8084C6D8 8FBF001C */  lw      $ra, 0x001C($sp)           
.L8084C6DC:
/* 1A4CC 8084C6DC C6040024 */  lwc1    $f4, 0x0024($s0)           ## 00000024
/* 1A4D0 8084C6E0 3C0141A0 */  lui     $at, 0x41A0                ## $at = 41A00000
/* 1A4D4 8084C6E4 44814000 */  mtc1    $at, $f8                   ## $f8 = 20.00
/* 1A4D8 8084C6E8 E7A40024 */  swc1    $f4, 0x0024($sp)           
/* 1A4DC 8084C6EC C6060028 */  lwc1    $f6, 0x0028($s0)           ## 00000028
/* 1A4E0 8084C6F0 8FA40044 */  lw      $a0, 0x0044($sp)           
/* 1A4E4 8084C6F4 27A50034 */  addiu   $a1, $sp, 0x0034           ## $a1 = FFFFFFF4
/* 1A4E8 8084C6F8 46083280 */  add.s   $f10, $f6, $f8             
/* 1A4EC 8084C6FC 248407C0 */  addiu   $a0, $a0, 0x07C0           ## $a0 = 000007C0
/* 1A4F0 8084C700 27A60030 */  addiu   $a2, $sp, 0x0030           ## $a2 = FFFFFFF0
/* 1A4F4 8084C704 27A70024 */  addiu   $a3, $sp, 0x0024           ## $a3 = FFFFFFE4
/* 1A4F8 8084C708 E7AA0028 */  swc1    $f10, 0x0028($sp)          
/* 1A4FC 8084C70C C610002C */  lwc1    $f16, 0x002C($s0)          ## 0000002C
/* 1A500 8084C710 AFA40020 */  sw      $a0, 0x0020($sp)           
/* 1A504 8084C714 0C00F250 */  jal     func_8003C940              
/* 1A508 8084C718 E7B0002C */  swc1    $f16, 0x002C($sp)          
/* 1A50C 8084C71C 44809000 */  mtc1    $zero, $f18                ## $f18 = 0.00
/* 1A510 8084C720 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 1A514 8084C724 8FA50034 */  lw      $a1, 0x0034($sp)           
/* 1A518 8084C728 46120032 */  c.eq.s  $f0, $f18                  
/* 1A51C 8084C72C 00000000 */  nop
/* 1A520 8084C730 45030007 */  bc1tl   .L8084C750                 
/* 1A524 8084C734 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 1A528 8084C738 0C0107C4 */  jal     func_80041F10              
/* 1A52C 8084C73C 8FA60030 */  lw      $a2, 0x0030($sp)           
/* 1A530 8084C740 A602089E */  sh      $v0, 0x089E($s0)           ## 0000089E
/* 1A534 8084C744 0C20CA28 */  jal     func_808328A0              
/* 1A538 8084C748 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
.L8084C74C:
/* 1A53C 8084C74C 8FBF001C */  lw      $ra, 0x001C($sp)           
.L8084C750:
/* 1A540 8084C750 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 1A544 8084C754 27BD0040 */  addiu   $sp, $sp, 0x0040           ## $sp = 00000000
/* 1A548 8084C758 03E00008 */  jr      $ra                        
/* 1A54C 8084C75C 00000000 */  nop


