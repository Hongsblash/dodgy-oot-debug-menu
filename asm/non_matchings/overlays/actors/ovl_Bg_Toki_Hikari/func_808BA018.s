glabel func_808BA018
/* 00118 808BA018 27BDFF98 */  addiu   $sp, $sp, 0xFF98           ## $sp = FFFFFF98
/* 0011C 808BA01C AFB10020 */  sw      $s1, 0x0020($sp)           
/* 00120 808BA020 00A08825 */  or      $s1, $a1, $zero            ## $s1 = 00000000
/* 00124 808BA024 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 00128 808BA028 AFB0001C */  sw      $s0, 0x001C($sp)           
/* 0012C 808BA02C AFA40068 */  sw      $a0, 0x0068($sp)           
/* 00130 808BA030 8CA50000 */  lw      $a1, 0x0000($a1)           ## 00000000
/* 00134 808BA034 3C06808C */  lui     $a2, %hi(D_808BAC80)       ## $a2 = 808C0000
/* 00138 808BA038 24C6AC80 */  addiu   $a2, $a2, %lo(D_808BAC80)  ## $a2 = 808BAC80
/* 0013C 808BA03C 27A40050 */  addiu   $a0, $sp, 0x0050           ## $a0 = FFFFFFE8
/* 00140 808BA040 240700F6 */  addiu   $a3, $zero, 0x00F6         ## $a3 = 000000F6
/* 00144 808BA044 0C031AB1 */  jal     func_800C6AC4              
/* 00148 808BA048 00A08025 */  or      $s0, $a1, $zero            ## $s0 = 00000000
/* 0014C 808BA04C 0C024F46 */  jal     func_80093D18              
/* 00150 808BA050 8E240000 */  lw      $a0, 0x0000($s1)           ## 00000000
/* 00154 808BA054 8E0202C0 */  lw      $v0, 0x02C0($s0)           ## 000002C0
/* 00158 808BA058 3C0FDA38 */  lui     $t7, 0xDA38                ## $t7 = DA380000
/* 0015C 808BA05C 35EF0003 */  ori     $t7, $t7, 0x0003           ## $t7 = DA380003
/* 00160 808BA060 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 00164 808BA064 AE0E02C0 */  sw      $t6, 0x02C0($s0)           ## 000002C0
/* 00168 808BA068 AC4F0000 */  sw      $t7, 0x0000($v0)           ## 00000000
/* 0016C 808BA06C 8E240000 */  lw      $a0, 0x0000($s1)           ## 00000000
/* 00170 808BA070 3C05808C */  lui     $a1, %hi(D_808BAC98)       ## $a1 = 808C0000
/* 00174 808BA074 24A5AC98 */  addiu   $a1, $a1, %lo(D_808BAC98)  ## $a1 = 808BAC98
/* 00178 808BA078 240600FC */  addiu   $a2, $zero, 0x00FC         ## $a2 = 000000FC
/* 0017C 808BA07C 0C0346A2 */  jal     Matrix_NewMtx              
/* 00180 808BA080 AFA2004C */  sw      $v0, 0x004C($sp)           
/* 00184 808BA084 8FA3004C */  lw      $v1, 0x004C($sp)           
/* 00188 808BA088 3C188016 */  lui     $t8, 0x8016                ## $t8 = 80160000
/* 0018C 808BA08C 3C09DE00 */  lui     $t1, 0xDE00                ## $t1 = DE000000
/* 00190 808BA090 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 00194 808BA094 8F18E664 */  lw      $t8, -0x199C($t8)          ## 8015E664
/* 00198 808BA098 3C0CDE00 */  lui     $t4, 0xDE00                ## $t4 = DE000000
/* 0019C 808BA09C 5700000A */  bnel    $t8, $zero, .L808BA0C8     
/* 001A0 808BA0A0 8E0202C0 */  lw      $v0, 0x02C0($s0)           ## 000002C0
/* 001A4 808BA0A4 8E0202C0 */  lw      $v0, 0x02C0($s0)           ## 000002C0
/* 001A8 808BA0A8 3C0A0601 */  lui     $t2, 0x0601                ## $t2 = 06010000
/* 001AC 808BA0AC 254A8190 */  addiu   $t2, $t2, 0x8190           ## $t2 = 06008190
/* 001B0 808BA0B0 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 001B4 808BA0B4 AE1902C0 */  sw      $t9, 0x02C0($s0)           ## 000002C0
/* 001B8 808BA0B8 AC4A0004 */  sw      $t2, 0x0004($v0)           ## 00000004
/* 001BC 808BA0BC 10000046 */  beq     $zero, $zero, .L808BA1D8   
/* 001C0 808BA0C0 AC490000 */  sw      $t1, 0x0000($v0)           ## 00000000
/* 001C4 808BA0C4 8E0202C0 */  lw      $v0, 0x02C0($s0)           ## 000002C0
.L808BA0C8:
/* 001C8 808BA0C8 3C0D0600 */  lui     $t5, 0x0600                ## $t5 = 06000000
/* 001CC 808BA0CC 25AD7E20 */  addiu   $t5, $t5, 0x7E20           ## $t5 = 06007E20
/* 001D0 808BA0D0 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 001D4 808BA0D4 AE0B02C0 */  sw      $t3, 0x02C0($s0)           ## 000002C0
/* 001D8 808BA0D8 AC4D0004 */  sw      $t5, 0x0004($v0)           ## 00000004
/* 001DC 808BA0DC AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
/* 001E0 808BA0E0 0C024F61 */  jal     func_80093D84              
/* 001E4 808BA0E4 8E240000 */  lw      $a0, 0x0000($s1)           ## 00000000
/* 001E8 808BA0E8 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 001EC 808BA0EC 3C0FDB06 */  lui     $t7, 0xDB06                ## $t7 = DB060000
/* 001F0 808BA0F0 35EF0020 */  ori     $t7, $t7, 0x0020           ## $t7 = DB060020
/* 001F4 808BA0F4 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 001F8 808BA0F8 AE0E02D0 */  sw      $t6, 0x02D0($s0)           ## 000002D0
/* 001FC 808BA0FC 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00200 808BA100 02214021 */  addu    $t0, $s1, $at              
/* 00204 808BA104 AC4F0000 */  sw      $t7, 0x0000($v0)           ## 00000000
/* 00208 808BA108 8D061DE4 */  lw      $a2, 0x1DE4($t0)           ## 00001DE4
/* 0020C 808BA10C 8E240000 */  lw      $a0, 0x0000($s1)           ## 00000000
/* 00210 808BA110 24180020 */  addiu   $t8, $zero, 0x0020         ## $t8 = 00000020
/* 00214 808BA114 AFB80010 */  sw      $t8, 0x0010($sp)           
/* 00218 808BA118 AFA80030 */  sw      $t0, 0x0030($sp)           
/* 0021C 808BA11C 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 00220 808BA120 24070040 */  addiu   $a3, $zero, 0x0040         ## $a3 = 00000040
/* 00224 808BA124 AFA20040 */  sw      $v0, 0x0040($sp)           
/* 00228 808BA128 0C0253A7 */  jal     Gfx_TexScroll              
/* 0022C 808BA12C 30C6007F */  andi    $a2, $a2, 0x007F           ## $a2 = 00000000
/* 00230 808BA130 8FA30040 */  lw      $v1, 0x0040($sp)           
/* 00234 808BA134 8FA80030 */  lw      $t0, 0x0030($sp)           
/* 00238 808BA138 3C09DB06 */  lui     $t1, 0xDB06                ## $t1 = DB060000
/* 0023C 808BA13C AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 00240 808BA140 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 00244 808BA144 35290024 */  ori     $t1, $t1, 0x0024           ## $t1 = DB060024
/* 00248 808BA148 240A0020 */  addiu   $t2, $zero, 0x0020         ## $t2 = 00000020
/* 0024C 808BA14C 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 00250 808BA150 AE1902D0 */  sw      $t9, 0x02D0($s0)           ## 000002D0
/* 00254 808BA154 AC490000 */  sw      $t1, 0x0000($v0)           ## 00000000
/* 00258 808BA158 8D061DE4 */  lw      $a2, 0x1DE4($t0)           ## 00001DE4
/* 0025C 808BA15C 8E240000 */  lw      $a0, 0x0000($s1)           ## 00000000
/* 00260 808BA160 AFAA0010 */  sw      $t2, 0x0010($sp)           
/* 00264 808BA164 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 00268 808BA168 24070040 */  addiu   $a3, $zero, 0x0040         ## $a3 = 00000040
/* 0026C 808BA16C AFA2003C */  sw      $v0, 0x003C($sp)           
/* 00270 808BA170 0C0253A7 */  jal     Gfx_TexScroll              
/* 00274 808BA174 30C6007F */  andi    $a2, $a2, 0x007F           ## $a2 = 00000000
/* 00278 808BA178 8FA3003C */  lw      $v1, 0x003C($sp)           
/* 0027C 808BA17C 3C0CDA38 */  lui     $t4, 0xDA38                ## $t4 = DA380000
/* 00280 808BA180 358C0003 */  ori     $t4, $t4, 0x0003           ## $t4 = DA380003
/* 00284 808BA184 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 00288 808BA188 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 0028C 808BA18C 3C05808C */  lui     $a1, %hi(D_808BACB0)       ## $a1 = 808C0000
/* 00290 808BA190 24A5ACB0 */  addiu   $a1, $a1, %lo(D_808BACB0)  ## $a1 = 808BACB0
/* 00294 808BA194 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 00298 808BA198 AE0B02D0 */  sw      $t3, 0x02D0($s0)           ## 000002D0
/* 0029C 808BA19C AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
/* 002A0 808BA1A0 8E240000 */  lw      $a0, 0x0000($s1)           ## 00000000
/* 002A4 808BA1A4 24060116 */  addiu   $a2, $zero, 0x0116         ## $a2 = 00000116
/* 002A8 808BA1A8 0C0346A2 */  jal     Matrix_NewMtx              
/* 002AC 808BA1AC AFA20038 */  sw      $v0, 0x0038($sp)           
/* 002B0 808BA1B0 8FA30038 */  lw      $v1, 0x0038($sp)           
/* 002B4 808BA1B4 3C0F0600 */  lui     $t7, 0x0600                ## $t7 = 06000000
/* 002B8 808BA1B8 25EF7EE0 */  addiu   $t7, $t7, 0x7EE0           ## $t7 = 06007EE0
/* 002BC 808BA1BC AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 002C0 808BA1C0 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 002C4 808BA1C4 3C0EDE00 */  lui     $t6, 0xDE00                ## $t6 = DE000000
/* 002C8 808BA1C8 244D0008 */  addiu   $t5, $v0, 0x0008           ## $t5 = 00000008
/* 002CC 808BA1CC AE0D02D0 */  sw      $t5, 0x02D0($s0)           ## 000002D0
/* 002D0 808BA1D0 AC4F0004 */  sw      $t7, 0x0004($v0)           ## 00000004
/* 002D4 808BA1D4 AC4E0000 */  sw      $t6, 0x0000($v0)           ## 00000000
.L808BA1D8:
/* 002D8 808BA1D8 3C06808C */  lui     $a2, %hi(D_808BACC8)       ## $a2 = 808C0000
/* 002DC 808BA1DC 24C6ACC8 */  addiu   $a2, $a2, %lo(D_808BACC8)  ## $a2 = 808BACC8
/* 002E0 808BA1E0 27A40050 */  addiu   $a0, $sp, 0x0050           ## $a0 = FFFFFFE8
/* 002E4 808BA1E4 8E250000 */  lw      $a1, 0x0000($s1)           ## 00000000
/* 002E8 808BA1E8 0C031AD5 */  jal     func_800C6B54              
/* 002EC 808BA1EC 2407011C */  addiu   $a3, $zero, 0x011C         ## $a3 = 0000011C
/* 002F0 808BA1F0 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 002F4 808BA1F4 8FB0001C */  lw      $s0, 0x001C($sp)           
/* 002F8 808BA1F8 8FB10020 */  lw      $s1, 0x0020($sp)           
/* 002FC 808BA1FC 03E00008 */  jr      $ra                        
/* 00300 808BA200 27BD0068 */  addiu   $sp, $sp, 0x0068           ## $sp = 00000000


