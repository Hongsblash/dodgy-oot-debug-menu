glabel func_80A791CC
/* 0021C 80A791CC AFA40000 */  sw      $a0, 0x0000($sp)           
/* 00220 80A791D0 94A2010E */  lhu     $v0, 0x010E($a1)           ## 0000010E
/* 00224 80A791D4 2401203E */  addiu   $at, $zero, 0x203E         ## $at = 0000203E
/* 00228 80A791D8 00001825 */  or      $v1, $zero, $zero          ## $v1 = 00000000
/* 0022C 80A791DC 1041000C */  beq     $v0, $at, .L80A79210       
/* 00230 80A791E0 2401203F */  addiu   $at, $zero, 0x203F         ## $at = 0000203F
/* 00234 80A791E4 1041000E */  beq     $v0, $at, .L80A79220       
/* 00238 80A791E8 24012045 */  addiu   $at, $zero, 0x2045         ## $at = 00002045
/* 0023C 80A791EC 14410014 */  bne     $v0, $at, .L80A79240       
/* 00240 80A791F0 3C028016 */  lui     $v0, 0x8016                ## $v0 = 80160000
/* 00244 80A791F4 2442E660 */  addiu   $v0, $v0, 0xE660           ## $v0 = 8015E660
/* 00248 80A791F8 944E0F0A */  lhu     $t6, 0x0F0A($v0)           ## 8015F56A
/* 0024C 80A791FC 35CF0080 */  ori     $t7, $t6, 0x0080           ## $t7 = 00000080
/* 00250 80A79200 A44F0F0A */  sh      $t7, 0x0F0A($v0)           ## 8015F56A
/* 00254 80A79204 00001400 */  sll     $v0, $zero, 16             
/* 00258 80A79208 03E00008 */  jr      $ra                        
/* 0025C 80A7920C 00021403 */  sra     $v0, $v0, 16               
.L80A79210:
/* 00260 80A79210 24030002 */  addiu   $v1, $zero, 0x0002         ## $v1 = 00000002
/* 00264 80A79214 00031400 */  sll     $v0, $v1, 16               
/* 00268 80A79218 03E00008 */  jr      $ra                        
/* 0026C 80A7921C 00021403 */  sra     $v0, $v0, 16               
.L80A79220:
/* 00270 80A79220 3C028016 */  lui     $v0, 0x8016                ## $v0 = 80160000
/* 00274 80A79224 2442E660 */  addiu   $v0, $v0, 0xE660           ## $v0 = 8015E660
/* 00278 80A79228 94580ED6 */  lhu     $t8, 0x0ED6($v0)           ## 8015F536
/* 0027C 80A7922C 94480F0A */  lhu     $t0, 0x0F0A($v0)           ## 8015F56A
/* 00280 80A79230 37190002 */  ori     $t9, $t8, 0x0002           ## $t9 = 00000002
/* 00284 80A79234 35090010 */  ori     $t1, $t0, 0x0010           ## $t1 = 00000010
/* 00288 80A79238 A4590ED6 */  sh      $t9, 0x0ED6($v0)           ## 8015F536
/* 0028C 80A7923C A4490F0A */  sh      $t1, 0x0F0A($v0)           ## 8015F56A
.L80A79240:
/* 00290 80A79240 00031400 */  sll     $v0, $v1, 16               
/* 00294 80A79244 03E00008 */  jr      $ra                        
/* 00298 80A79248 00021403 */  sra     $v0, $v0, 16               


