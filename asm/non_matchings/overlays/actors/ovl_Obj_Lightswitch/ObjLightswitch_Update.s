glabel ObjLightswitch_Update
/* 00B4C 80B976EC 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 00B50 80B976F0 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 00B54 80B976F4 AFB10018 */  sw      $s1, 0x0018($sp)           
/* 00B58 80B976F8 AFB00014 */  sw      $s0, 0x0014($sp)           
/* 00B5C 80B976FC 848201B2 */  lh      $v0, 0x01B2($a0)           ## 000001B2
/* 00B60 80B97700 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00B64 80B97704 00A08825 */  or      $s1, $a1, $zero            ## $s1 = 00000000
/* 00B68 80B97708 18400002 */  blez    $v0, .L80B97714            
/* 00B6C 80B9770C 244EFFFF */  addiu   $t6, $v0, 0xFFFF           ## $t6 = FFFFFFFF
/* 00B70 80B97710 A48E01B2 */  sh      $t6, 0x01B2($a0)           ## 000001B2
.L80B97714:
/* 00B74 80B97714 8E19014C */  lw      $t9, 0x014C($s0)           ## 0000014C
/* 00B78 80B97718 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00B7C 80B9771C 02202825 */  or      $a1, $s1, $zero            ## $a1 = 00000000
/* 00B80 80B97720 0320F809 */  jalr    $ra, $t9                   
/* 00B84 80B97724 00000000 */  nop
/* 00B88 80B97728 8E0F0130 */  lw      $t7, 0x0130($s0)           ## 00000130
/* 00B8C 80B9772C 51E00024 */  beql    $t7, $zero, .L80B977C0     
/* 00B90 80B97730 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 00B94 80B97734 8618001C */  lh      $t8, 0x001C($s0)           ## 0000001C
/* 00B98 80B97738 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 00B9C 80B9773C 33080001 */  andi    $t0, $t8, 0x0001           ## $t0 = 00000000
/* 00BA0 80B97740 5501000F */  bnel    $t0, $at, .L80B97780       
/* 00BA4 80B97744 92020161 */  lbu     $v0, 0x0161($s0)           ## 00000161
/* 00BA8 80B97748 8E02011C */  lw      $v0, 0x011C($s0)           ## 0000011C
/* 00BAC 80B9774C 3C014270 */  lui     $at, 0x4270                ## $at = 42700000
/* 00BB0 80B97750 44814000 */  mtc1    $at, $f8                   ## $f8 = 60.00
/* 00BB4 80B97754 C4440024 */  lwc1    $f4, 0x0024($v0)           ## 00000024
/* 00BB8 80B97758 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00BBC 80B9775C 24050000 */  addiu   $a1, $zero, 0x0000         ## $a1 = 00000000
/* 00BC0 80B97760 E6040024 */  swc1    $f4, 0x0024($s0)           ## 00000024
/* 00BC4 80B97764 C4460028 */  lwc1    $f6, 0x0028($v0)           ## 00000028
/* 00BC8 80B97768 46083280 */  add.s   $f10, $f6, $f8             
/* 00BCC 80B9776C E60A0028 */  swc1    $f10, 0x0028($s0)          ## 00000028
/* 00BD0 80B97770 C450002C */  lwc1    $f16, 0x002C($v0)          ## 0000002C
/* 00BD4 80B97774 0C00B56E */  jal     Actor_SetHeight
              
/* 00BD8 80B97778 E610002C */  swc1    $f16, 0x002C($s0)          ## 0000002C
/* 00BDC 80B9777C 92020161 */  lbu     $v0, 0x0161($s0)           ## 00000161
.L80B97780:
/* 00BE0 80B97780 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00BE4 80B97784 34211E60 */  ori     $at, $at, 0x1E60           ## $at = 00011E60
/* 00BE8 80B97788 3049FFFD */  andi    $t1, $v0, 0xFFFD           ## $t1 = 00000000
/* 00BEC 80B9778C A2090161 */  sb      $t1, 0x0161($s0)           ## 00000161
/* 00BF0 80B97790 02212821 */  addu    $a1, $s1, $at              
/* 00BF4 80B97794 26060150 */  addiu   $a2, $s0, 0x0150           ## $a2 = 00000150
/* 00BF8 80B97798 A20201C2 */  sb      $v0, 0x01C2($s0)           ## 000001C2
/* 00BFC 80B9779C AFA60020 */  sw      $a2, 0x0020($sp)           
/* 00C00 80B977A0 AFA50024 */  sw      $a1, 0x0024($sp)           
/* 00C04 80B977A4 0C017713 */  jal     Actor_CollisionCheck_SetOT
              ## CollisionCheck_setOT
/* 00C08 80B977A8 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 00C0C 80B977AC 8FA50024 */  lw      $a1, 0x0024($sp)           
/* 00C10 80B977B0 8FA60020 */  lw      $a2, 0x0020($sp)           
/* 00C14 80B977B4 0C01767D */  jal     Actor_CollisionCheck_SetAC
              ## CollisionCheck_setAC
/* 00C18 80B977B8 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 00C1C 80B977BC 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80B977C0:
/* 00C20 80B977C0 8FB00014 */  lw      $s0, 0x0014($sp)           
/* 00C24 80B977C4 8FB10018 */  lw      $s1, 0x0018($sp)           
/* 00C28 80B977C8 03E00008 */  jr      $ra                        
/* 00C2C 80B977CC 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000


