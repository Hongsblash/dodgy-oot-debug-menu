glabel func_80B93BF0
/* 00090 80B93BF0 27BDFFC0 */  addiu   $sp, $sp, 0xFFC0           ## $sp = FFFFFFC0
/* 00094 80B93BF4 AFBF003C */  sw      $ra, 0x003C($sp)           
/* 00098 80B93BF8 AFB00038 */  sw      $s0, 0x0038($sp)           
/* 0009C 80B93BFC 848E001C */  lh      $t6, 0x001C($a0)           ## 0000001C
/* 000A0 80B93C00 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 000A4 80B93C04 00A03025 */  or      $a2, $a1, $zero            ## $a2 = 00000000
/* 000A8 80B93C08 000E7943 */  sra     $t7, $t6,  5               
/* 000AC 80B93C0C 31F80001 */  andi    $t8, $t7, 0x0001           ## $t8 = 00000000
/* 000B0 80B93C10 53000014 */  beql    $t8, $zero, .L80B93C64     
/* 000B4 80B93C14 8FBF003C */  lw      $ra, 0x003C($sp)           
/* 000B8 80B93C18 C6040024 */  lwc1    $f4, 0x0024($s0)           ## 00000024
/* 000BC 80B93C1C 24A41C24 */  addiu   $a0, $a1, 0x1C24           ## $a0 = 00001C24
/* 000C0 80B93C20 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 000C4 80B93C24 E7A40010 */  swc1    $f4, 0x0010($sp)           
/* 000C8 80B93C28 C6060028 */  lwc1    $f6, 0x0028($s0)           ## 00000028
/* 000CC 80B93C2C 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 000D0 80B93C30 2407011E */  addiu   $a3, $zero, 0x011E         ## $a3 = 0000011E
/* 000D4 80B93C34 E7A60014 */  swc1    $f6, 0x0014($sp)           
/* 000D8 80B93C38 C608002C */  lwc1    $f8, 0x002C($s0)           ## 0000002C
/* 000DC 80B93C3C E7A80018 */  swc1    $f8, 0x0018($sp)           
/* 000E0 80B93C40 86190030 */  lh      $t9, 0x0030($s0)           ## 00000030
/* 000E4 80B93C44 AFB9001C */  sw      $t9, 0x001C($sp)           
/* 000E8 80B93C48 86080032 */  lh      $t0, 0x0032($s0)           ## 00000032
/* 000EC 80B93C4C AFA80020 */  sw      $t0, 0x0020($sp)           
/* 000F0 80B93C50 86090034 */  lh      $t1, 0x0034($s0)           ## 00000034
/* 000F4 80B93C54 AFAA0028 */  sw      $t2, 0x0028($sp)           
/* 000F8 80B93C58 0C00C916 */  jal     Actor_SpawnAttached
              
/* 000FC 80B93C5C AFA90024 */  sw      $t1, 0x0024($sp)           
/* 00100 80B93C60 8FBF003C */  lw      $ra, 0x003C($sp)           
.L80B93C64:
/* 00104 80B93C64 8FB00038 */  lw      $s0, 0x0038($sp)           
/* 00108 80B93C68 27BD0040 */  addiu   $sp, $sp, 0x0040           ## $sp = 00000000
/* 0010C 80B93C6C 03E00008 */  jr      $ra                        
/* 00110 80B93C70 00000000 */  nop


