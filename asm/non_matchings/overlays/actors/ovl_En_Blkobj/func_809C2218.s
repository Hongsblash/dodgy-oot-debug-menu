glabel func_809C2218
/* 001B8 809C2218 27BDFFE0 */  addiu   $sp, $sp, 0xFFE0           ## $sp = FFFFFFE0
/* 001BC 809C221C AFBF001C */  sw      $ra, 0x001C($sp)           
/* 001C0 809C2220 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 001C4 809C2224 84820166 */  lh      $v0, 0x0166($a0)           ## 00000166
/* 001C8 809C2228 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 001CC 809C222C 00A03825 */  or      $a3, $a1, $zero            ## $a3 = 00000000
/* 001D0 809C2230 1440000F */  bne     $v0, $zero, .L809C2270     
/* 001D4 809C2234 28430065 */  slti    $v1, $v0, 0x0065           
/* 001D8 809C2238 24A41C24 */  addiu   $a0, $a1, 0x1C24           ## $a0 = 00001C24
/* 001DC 809C223C 24050033 */  addiu   $a1, $zero, 0x0033         ## $a1 = 00000033
/* 001E0 809C2240 24060009 */  addiu   $a2, $zero, 0x0009         ## $a2 = 00000009
/* 001E4 809C2244 0C00CB0F */  jal     Actor_Find
              
/* 001E8 809C2248 AFA70024 */  sw      $a3, 0x0024($sp)           
/* 001EC 809C224C 14400024 */  bne     $v0, $zero, .L809C22E0     
/* 001F0 809C2250 8FA70024 */  lw      $a3, 0x0024($sp)           
/* 001F4 809C2254 00E02025 */  or      $a0, $a3, $zero            ## $a0 = 00000000
/* 001F8 809C2258 0C00B33C */  jal     Flags_SetClear
              
/* 001FC 809C225C 82050003 */  lb      $a1, 0x0003($s0)           ## 00000003
/* 00200 809C2260 860E0166 */  lh      $t6, 0x0166($s0)           ## 00000166
/* 00204 809C2264 25CF0001 */  addiu   $t7, $t6, 0x0001           ## $t7 = 00000001
/* 00208 809C2268 1000001D */  beq     $zero, $zero, .L809C22E0   
/* 0020C 809C226C A60F0166 */  sh      $t7, 0x0166($s0)           ## 00000166
.L809C2270:
/* 00210 809C2270 38630001 */  xori    $v1, $v1, 0x0001           ## $v1 = 00000001
/* 00214 809C2274 24580001 */  addiu   $t8, $v0, 0x0001           ## $t8 = 00000001
/* 00218 809C2278 10600019 */  beq     $v1, $zero, .L809C22E0     
/* 0021C 809C227C A6180166 */  sh      $t8, 0x0166($s0)           ## 00000166
/* 00220 809C2280 86020166 */  lh      $v0, 0x0166($s0)           ## 00000166
/* 00224 809C2284 240A00FF */  addiu   $t2, $zero, 0x00FF         ## $t2 = 000000FF
/* 00228 809C2288 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 0022C 809C228C 2442FF9C */  addiu   $v0, $v0, 0xFF9C           ## $v0 = FFFFFF9C
/* 00230 809C2290 00021083 */  sra     $v0, $v0,  2               
/* 00234 809C2294 28410006 */  slti    $at, $v0, 0x0006           
/* 00238 809C2298 14200002 */  bne     $at, $zero, .L809C22A4     
/* 0023C 809C229C 3C05809C */  lui     $a1, %hi(func_809C22F4)    ## $a1 = 809C0000
/* 00240 809C22A0 24020005 */  addiu   $v0, $zero, 0x0005         ## $v0 = 00000005
.L809C22A4:
/* 00244 809C22A4 86190164 */  lh      $t9, 0x0164($s0)           ## 00000164
/* 00248 809C22A8 24A522F4 */  addiu   $a1, $a1, %lo(func_809C22F4) ## $a1 = 809C22F4
/* 0024C 809C22AC 03224021 */  addu    $t0, $t9, $v0              
/* 00250 809C22B0 A6080164 */  sh      $t0, 0x0164($s0)           ## 00000164
/* 00254 809C22B4 86090164 */  lh      $t1, 0x0164($s0)           ## 00000164
/* 00258 809C22B8 29210100 */  slti    $at, $t1, 0x0100           
/* 0025C 809C22BC 54200009 */  bnel    $at, $zero, .L809C22E4     
/* 00260 809C22C0 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 00264 809C22C4 A60A0164 */  sh      $t2, 0x0164($s0)           ## 00000164
/* 00268 809C22C8 0C270818 */  jal     func_809C2060              
/* 0026C 809C22CC AFA70024 */  sw      $a3, 0x0024($sp)           
/* 00270 809C22D0 8FA40024 */  lw      $a0, 0x0024($sp)           
/* 00274 809C22D4 8E06014C */  lw      $a2, 0x014C($s0)           ## 0000014C
/* 00278 809C22D8 0C00FB56 */  jal     DynaPolyInfo_Free
              ## DynaPolyInfo_delReserve
/* 0027C 809C22DC 24850810 */  addiu   $a1, $a0, 0x0810           ## $a1 = 00000810
.L809C22E0:
/* 00280 809C22E0 8FBF001C */  lw      $ra, 0x001C($sp)           
.L809C22E4:
/* 00284 809C22E4 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 00288 809C22E8 27BD0020 */  addiu   $sp, $sp, 0x0020           ## $sp = 00000000
/* 0028C 809C22EC 03E00008 */  jr      $ra                        
/* 00290 809C22F0 00000000 */  nop


