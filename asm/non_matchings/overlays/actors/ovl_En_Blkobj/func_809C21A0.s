glabel func_809C21A0
/* 00140 809C21A0 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 00144 809C21A4 AFBF0034 */  sw      $ra, 0x0034($sp)           
/* 00148 809C21A8 AFB00030 */  sw      $s0, 0x0030($sp)           
/* 0014C 809C21AC 8C8E0004 */  lw      $t6, 0x0004($a0)           ## 00000004
/* 00150 809C21B0 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00154 809C21B4 31CF0040 */  andi    $t7, $t6, 0x0040           ## $t7 = 00000000
/* 00158 809C21B8 55E00013 */  bnel    $t7, $zero, .L809C2208     
/* 0015C 809C21BC 8FBF0034 */  lw      $ra, 0x0034($sp)           
/* 00160 809C21C0 C6040028 */  lwc1    $f4, 0x0028($s0)           ## 00000028
/* 00164 809C21C4 8E070024 */  lw      $a3, 0x0024($s0)           ## 00000024
/* 00168 809C21C8 24A41C24 */  addiu   $a0, $a1, 0x1C24           ## $a0 = 00001C24
/* 0016C 809C21CC E7A40010 */  swc1    $f4, 0x0010($sp)           
/* 00170 809C21D0 C606002C */  lwc1    $f6, 0x002C($s0)           ## 0000002C
/* 00174 809C21D4 AFA00018 */  sw      $zero, 0x0018($sp)         
/* 00178 809C21D8 24060033 */  addiu   $a2, $zero, 0x0033         ## $a2 = 00000033
/* 0017C 809C21DC E7A60014 */  swc1    $f6, 0x0014($sp)           
/* 00180 809C21E0 8618008A */  lh      $t8, 0x008A($s0)           ## 0000008A
/* 00184 809C21E4 AFA00024 */  sw      $zero, 0x0024($sp)         
/* 00188 809C21E8 AFA00020 */  sw      $zero, 0x0020($sp)         
/* 0018C 809C21EC 0C00C7D4 */  jal     Actor_Spawn
              ## ActorSpawn
/* 00190 809C21F0 AFB8001C */  sw      $t8, 0x001C($sp)           
/* 00194 809C21F4 3C05809C */  lui     $a1, %hi(func_809C2218)    ## $a1 = 809C0000
/* 00198 809C21F8 24A52218 */  addiu   $a1, $a1, %lo(func_809C2218) ## $a1 = 809C2218
/* 0019C 809C21FC 0C270818 */  jal     func_809C2060              
/* 001A0 809C2200 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 001A4 809C2204 8FBF0034 */  lw      $ra, 0x0034($sp)           
.L809C2208:
/* 001A8 809C2208 8FB00030 */  lw      $s0, 0x0030($sp)           
/* 001AC 809C220C 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 001B0 809C2210 03E00008 */  jr      $ra                        
/* 001B4 809C2214 00000000 */  nop
