glabel func_80A710F8
/* 01B48 80A710F8 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 01B4C 80A710FC AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 01B50 80A71100 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 01B54 80A71104 848E01E8 */  lh      $t6, 0x01E8($a0)           ## 000001E8
/* 01B58 80A71108 00803825 */  or      $a3, $a0, $zero            ## $a3 = 00000000
/* 01B5C 80A7110C 3C198016 */  lui     $t9, %hi(gSaveContext+0x1400)
/* 01B60 80A71110 11C0000C */  beq     $t6, $zero, .L80A71144     
/* 01B64 80A71114 00000000 */  nop
/* 01B68 80A71118 8C980154 */  lw      $t8, 0x0154($a0)           ## 00000154
/* 01B6C 80A7111C 3C0F0600 */  lui     $t7, 0x0600                ## $t7 = 06000000
/* 01B70 80A71120 25EF0BFC */  addiu   $t7, $t7, 0x0BFC           ## $t7 = 06000BFC
/* 01B74 80A71124 11F8001F */  beq     $t7, $t8, .L80A711A4       
/* 01B78 80A71128 2484014C */  addiu   $a0, $a0, 0x014C           ## $a0 = 0000014C
/* 01B7C 80A7112C 3C0580A7 */  lui     $a1, %hi(D_80A72050)       ## $a1 = 80A70000
/* 01B80 80A71130 24A52050 */  addiu   $a1, $a1, %lo(D_80A72050)  ## $a1 = 80A72050
/* 01B84 80A71134 0C00D3B0 */  jal     func_80034EC0              
/* 01B88 80A71138 2406001A */  addiu   $a2, $zero, 0x001A         ## $a2 = 0000001A
/* 01B8C 80A7113C 1000001A */  beq     $zero, $zero, .L80A711A8   
/* 01B90 80A71140 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80A71144:
/* 01B94 80A71144 9739FA60 */  lhu     $t9, %lo(gSaveContext+0x1400)($t9)
/* 01B98 80A71148 33280001 */  andi    $t0, $t9, 0x0001           ## $t0 = 00000000
/* 01B9C 80A7114C 5100000D */  beql    $t0, $zero, .L80A71184     
/* 01BA0 80A71150 8CEC0154 */  lw      $t4, 0x0154($a3)           ## 00000154
/* 01BA4 80A71154 8CEA0154 */  lw      $t2, 0x0154($a3)           ## 00000154
/* 01BA8 80A71158 3C090600 */  lui     $t1, 0x0600                ## $t1 = 06000000
/* 01BAC 80A7115C 25290FE4 */  addiu   $t1, $t1, 0x0FE4           ## $t1 = 06000FE4
/* 01BB0 80A71160 112A0010 */  beq     $t1, $t2, .L80A711A4       
/* 01BB4 80A71164 24E4014C */  addiu   $a0, $a3, 0x014C           ## $a0 = 0000014C
/* 01BB8 80A71168 3C0580A7 */  lui     $a1, %hi(D_80A72050)       ## $a1 = 80A70000
/* 01BBC 80A7116C 24A52050 */  addiu   $a1, $a1, %lo(D_80A72050)  ## $a1 = 80A72050
/* 01BC0 80A71170 0C00D3B0 */  jal     func_80034EC0              
/* 01BC4 80A71174 24060019 */  addiu   $a2, $zero, 0x0019         ## $a2 = 00000019
/* 01BC8 80A71178 1000000B */  beq     $zero, $zero, .L80A711A8   
/* 01BCC 80A7117C 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 01BD0 80A71180 8CEC0154 */  lw      $t4, 0x0154($a3)           ## 00000154
.L80A71184:
/* 01BD4 80A71184 3C0B0600 */  lui     $t3, 0x0600                ## $t3 = 06000000
/* 01BD8 80A71188 256B12E8 */  addiu   $t3, $t3, 0x12E8           ## $t3 = 060012E8
/* 01BDC 80A7118C 116C0005 */  beq     $t3, $t4, .L80A711A4       
/* 01BE0 80A71190 24E4014C */  addiu   $a0, $a3, 0x014C           ## $a0 = 0000014C
/* 01BE4 80A71194 3C0580A7 */  lui     $a1, %hi(D_80A72050)       ## $a1 = 80A70000
/* 01BE8 80A71198 24A52050 */  addiu   $a1, $a1, %lo(D_80A72050)  ## $a1 = 80A72050
/* 01BEC 80A7119C 0C00D3B0 */  jal     func_80034EC0              
/* 01BF0 80A711A0 24060018 */  addiu   $a2, $zero, 0x0018         ## $a2 = 00000018
.L80A711A4:
/* 01BF4 80A711A4 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80A711A8:
/* 01BF8 80A711A8 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 01BFC 80A711AC 03E00008 */  jr      $ra                        
/* 01C00 80A711B0 00000000 */  nop
