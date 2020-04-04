glabel BgIngate_Draw
/* 00280 808929C0 27BDFFC0 */  addiu   $sp, $sp, 0xFFC0           ## $sp = FFFFFFC0
/* 00284 808929C4 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 00288 808929C8 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 0028C 808929CC AFA40040 */  sw      $a0, 0x0040($sp)           
/* 00290 808929D0 AFA50044 */  sw      $a1, 0x0044($sp)           
/* 00294 808929D4 8CA50000 */  lw      $a1, 0x0000($a1)           ## 00000000
/* 00298 808929D8 3C068089 */  lui     $a2, %hi(D_80892AB0)       ## $a2 = 80890000
/* 0029C 808929DC 24C62AB0 */  addiu   $a2, $a2, %lo(D_80892AB0)  ## $a2 = 80892AB0
/* 002A0 808929E0 27A4002C */  addiu   $a0, $sp, 0x002C           ## $a0 = FFFFFFEC
/* 002A4 808929E4 240700F0 */  addiu   $a3, $zero, 0x00F0         ## $a3 = 000000F0
/* 002A8 808929E8 0C031AB1 */  jal     Graph_OpenDisp              
/* 002AC 808929EC 00A08025 */  or      $s0, $a1, $zero            ## $s0 = 00000000
/* 002B0 808929F0 8FAF0044 */  lw      $t7, 0x0044($sp)           
/* 002B4 808929F4 0C024F46 */  jal     func_80093D18              
/* 002B8 808929F8 8DE40000 */  lw      $a0, 0x0000($t7)           ## 00000000
/* 002BC 808929FC 8E0202C0 */  lw      $v0, 0x02C0($s0)           ## 000002C0
/* 002C0 80892A00 3C19DA38 */  lui     $t9, 0xDA38                ## $t9 = DA380000
/* 002C4 80892A04 37390003 */  ori     $t9, $t9, 0x0003           ## $t9 = DA380003
/* 002C8 80892A08 24580008 */  addiu   $t8, $v0, 0x0008           ## $t8 = 00000008
/* 002CC 80892A0C AE1802C0 */  sw      $t8, 0x02C0($s0)           ## 000002C0
/* 002D0 80892A10 AC590000 */  sw      $t9, 0x0000($v0)           ## 00000000
/* 002D4 80892A14 8FA80044 */  lw      $t0, 0x0044($sp)           
/* 002D8 80892A18 3C058089 */  lui     $a1, %hi(D_80892AC4)       ## $a1 = 80890000
/* 002DC 80892A1C 24A52AC4 */  addiu   $a1, $a1, %lo(D_80892AC4)  ## $a1 = 80892AC4
/* 002E0 80892A20 8D040000 */  lw      $a0, 0x0000($t0)           ## 00000000
/* 002E4 80892A24 240600F5 */  addiu   $a2, $zero, 0x00F5         ## $a2 = 000000F5
/* 002E8 80892A28 0C0346A2 */  jal     Matrix_NewMtx              
/* 002EC 80892A2C AFA20028 */  sw      $v0, 0x0028($sp)           
/* 002F0 80892A30 8FA30028 */  lw      $v1, 0x0028($sp)           
/* 002F4 80892A34 3C0B0600 */  lui     $t3, 0x0600                ## $t3 = 06000000
/* 002F8 80892A38 256B1040 */  addiu   $t3, $t3, 0x1040           ## $t3 = 06001040
/* 002FC 80892A3C AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 00300 80892A40 8E0202C0 */  lw      $v0, 0x02C0($s0)           ## 000002C0
/* 00304 80892A44 3C0ADE00 */  lui     $t2, 0xDE00                ## $t2 = DE000000
/* 00308 80892A48 3C068089 */  lui     $a2, %hi(D_80892AD8)       ## $a2 = 80890000
/* 0030C 80892A4C 24490008 */  addiu   $t1, $v0, 0x0008           ## $t1 = 00000008
/* 00310 80892A50 AE0902C0 */  sw      $t1, 0x02C0($s0)           ## 000002C0
/* 00314 80892A54 AC4B0004 */  sw      $t3, 0x0004($v0)           ## 00000004
/* 00318 80892A58 AC4A0000 */  sw      $t2, 0x0000($v0)           ## 00000000
/* 0031C 80892A5C 8FAC0044 */  lw      $t4, 0x0044($sp)           
/* 00320 80892A60 24C62AD8 */  addiu   $a2, $a2, %lo(D_80892AD8)  ## $a2 = 80892AD8
/* 00324 80892A64 27A4002C */  addiu   $a0, $sp, 0x002C           ## $a0 = FFFFFFEC
/* 00328 80892A68 240700FA */  addiu   $a3, $zero, 0x00FA         ## $a3 = 000000FA
/* 0032C 80892A6C 0C031AD5 */  jal     Graph_CloseDisp              
/* 00330 80892A70 8D850000 */  lw      $a1, 0x0000($t4)           ## 00000000
/* 00334 80892A74 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 00338 80892A78 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 0033C 80892A7C 27BD0040 */  addiu   $sp, $sp, 0x0040           ## $sp = 00000000
/* 00340 80892A80 03E00008 */  jr      $ra                        
/* 00344 80892A84 00000000 */  nop
/* 00348 80892A88 00000000 */  nop
/* 0034C 80892A8C 00000000 */  nop

