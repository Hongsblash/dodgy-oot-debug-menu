glabel func_80846050
/* 13E40 80846050 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 13E44 80846054 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 13E48 80846058 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 13E4C 8084605C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 13E50 80846060 0C20DC87 */  jal     func_8083721C              
/* 13E54 80846064 AFA5002C */  sw      $a1, 0x002C($sp)           
/* 13E58 80846068 260601B4 */  addiu   $a2, $s0, 0x01B4           ## $a2 = 000001B4
/* 13E5C 8084606C 00C02825 */  or      $a1, $a2, $zero            ## $a1 = 000001B4
/* 13E60 80846070 AFA60020 */  sw      $a2, 0x0020($sp)           
/* 13E64 80846074 0C028EF0 */  jal     func_800A3BC0              
/* 13E68 80846078 8FA4002C */  lw      $a0, 0x002C($sp)           
/* 13E6C 8084607C 10400009 */  beq     $v0, $zero, .L808460A4     
/* 13E70 80846080 8FA60020 */  lw      $a2, 0x0020($sp)           
/* 13E74 80846084 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 13E78 80846088 0C20E7E4 */  jal     func_80839F90              
/* 13E7C 8084608C 8FA5002C */  lw      $a1, 0x002C($sp)           
/* 13E80 80846090 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 13E84 80846094 0C20D5A2 */  jal     func_80835688              
/* 13E88 80846098 8FA5002C */  lw      $a1, 0x002C($sp)           
/* 13E8C 8084609C 1000001C */  beq     $zero, $zero, .L80846110   
/* 13E90 808460A0 8FBF001C */  lw      $ra, 0x001C($sp)           
.L808460A4:
/* 13E94 808460A4 00C02025 */  or      $a0, $a2, $zero            ## $a0 = 00000000
/* 13E98 808460A8 0C02914C */  jal     func_800A4530              
/* 13E9C 808460AC 3C054080 */  lui     $a1, 0x4080                ## $a1 = 40800000
/* 13EA0 808460B0 10400013 */  beq     $v0, $zero, .L80846100     
/* 13EA4 808460B4 260403BE */  addiu   $a0, $s0, 0x03BE           ## $a0 = 000003BE
/* 13EA8 808460B8 8E060438 */  lw      $a2, 0x0438($s0)           ## 00000438
/* 13EAC 808460BC 8FA4002C */  lw      $a0, 0x002C($sp)           
/* 13EB0 808460C0 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 13EB4 808460C4 0C20D591 */  jal     func_80835644              
/* 13EB8 808460C8 AFA60024 */  sw      $a2, 0x0024($sp)           
/* 13EBC 808460CC 1440000F */  bne     $v0, $zero, .L8084610C     
/* 13EC0 808460D0 8FA60024 */  lw      $a2, 0x0024($sp)           
/* 13EC4 808460D4 AE0603AC */  sw      $a2, 0x03AC($s0)           ## 000003AC
/* 13EC8 808460D8 AE06011C */  sw      $a2, 0x011C($s0)           ## 0000011C
/* 13ECC 808460DC 94CE0088 */  lhu     $t6, 0x0088($a2)           ## 00000088
/* 13ED0 808460E0 ACD00118 */  sw      $s0, 0x0118($a2)           ## 00000118
/* 13ED4 808460E4 84D800B6 */  lh      $t8, 0x00B6($a2)           ## 000000B6
/* 13ED8 808460E8 31CFFF00 */  andi    $t7, $t6, 0xFF00           ## $t7 = 00000000
/* 13EDC 808460EC A4CF0088 */  sh      $t7, 0x0088($a2)           ## 00000088
/* 13EE0 808460F0 861900B6 */  lh      $t9, 0x00B6($s0)           ## 000000B6
/* 13EE4 808460F4 03194023 */  subu    $t0, $t8, $t9              
/* 13EE8 808460F8 10000004 */  beq     $zero, $zero, .L8084610C   
/* 13EEC 808460FC A60803BE */  sh      $t0, 0x03BE($s0)           ## 000003BE
.L80846100:
/* 13EF0 80846100 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 13EF4 80846104 0C01DE2B */  jal     Math_ApproxUpdateScaledS
              
/* 13EF8 80846108 24060FA0 */  addiu   $a2, $zero, 0x0FA0         ## $a2 = 00000FA0
.L8084610C:
/* 13EFC 8084610C 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80846110:
/* 13F00 80846110 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 13F04 80846114 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 13F08 80846118 03E00008 */  jr      $ra                        
/* 13F0C 8084611C 00000000 */  nop


