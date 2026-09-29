'use strict';
// 手がかりを京大の問題文の書き方（観察の文）にし、その結果から「どの部分構造が決まるか」を言葉にする。
// sentence: 問題文に載せる文、infer: ヒント（この実験からわかること）

// 観察の文。subject は「化合物 A」など。result はカードの結果
// products のときは文の最後に「次の化合物が得られた」と書き、構造式は画面側で並べる
const S = {
  silver_mirror: (s, r) => `${s} をアンモニア性硝酸銀水溶液と温めると、${r ? '銀が析出した' : '銀は析出しなかった'}。`,
  fehling: (s, r) => `${s} をフェーリング液と加熱すると、${r ? '赤色の沈殿が生じた' : '変化はなかった'}。`,
  iodoform: (s, r) => `${s} に水酸化ナトリウム水溶液とヨウ素を加えて温めると、${r ? '特有のにおいをもつ黄色沈殿が生じた' : '沈殿は生じなかった'}。`,
  sodium: (s, r) => `${s} に金属ナトリウムを加えると、${r ? '気体が発生した' : '変化はなかった'}。`,
  nahco3: (s, r) => `${s} は炭酸水素ナトリウム水溶液に${r ? '気体を発生しながら溶けた' : '溶けなかった'}。`,
  fecl3: (s, r) => `${s} に塩化鉄(III) 水溶液を加えると、${r ? '呈色した' : '呈色しなかった'}。`,
  bromine: (s, r) => `${s} は臭素水の色を${r ? '消した' : '消さなかった'}。`,
  naoh: (s, r) => `${s} は水酸化ナトリウム水溶液に${r ? '溶けた' : '溶けなかった'}。`,
  hcl: (s, r) => `${s} は希塩酸に${r ? '溶けた' : '溶けなかった'}。`,
  chiral: (s, r) => (r ? `${s} は不斉炭素原子を ${r} 個もつ。` : `${s} は不斉炭素原子をもたない。`),
  cis_trans: (s, r) => `${s} にはシス-トランス異性体が${r ? '存在する' : '存在しない'}。`,
  anhydride: (s, r) => `${s} を加熱すると、${r ? '分子内で脱水して酸無水物が生じた' : '酸無水物は生じなかった'}。`,
  carbon_env: (s, r) => `${s} には、化学的に等価でない炭素原子が ${r} 種類ある。`,
  ring_cl: (s, r) => `${s} のベンゼン環の水素原子1個を塩素原子で置き換えた化合物は、${r} 種類考えられる。`,
  cl_sub: (s, r) => `${s} の炭素原子に結合した水素原子1個を塩素原子で置き換えた化合物は、構造異性体として ${r} 種類考えられる。`,
  h2_uptake: (s, r) => `${s} 1 mol に白金触媒を用いて水素を付加させると、${r} mol の水素が消費された。`,
  dehydration_count: (s, r) => `${s} を脱水して得られるアルケンは、シス-トランス異性体を区別して ${r} 種類であった。`,
  stereo_count: (s, r) => `${s} には、${s.replace(/^化合物 /, '')} 自身を含めて ${r} 種類の立体異性体が存在する。`,
  optically_active: (s, r) => `${s} の水溶液は偏光面を${r ? '回転させた' : '回転させなかった'}。`,
  ninhydrin: (s, r) => `${s} にニンヒドリン水溶液を加えて温めると、${r ? '赤紫色になった' : '変化はなかった'}。`,
  xanthoprotein: (s, r) => `${s} に濃硝酸を加えて加熱すると、${r ? '黄色になった' : '変化はなかった'}。`,
  sulfur: (s, r) => `${s} に水酸化ナトリウムを加えて加熱し、酢酸鉛(II) 水溶液を加えると、${r ? '黒色沈殿が生じた' : '変化はなかった'}。`,
  biuret: (s, r) => `${s} にビウレット反応を行うと、${r ? '赤紫色になった' : '青色のままだった'}。`,
  alpha_amino: (s, r) => `${s} は α-アミノ酸${r ? 'である' : 'ではない'}。`,
};
// 生成物を返す反応の言い回し（「〜すると」の部分）
const V = {
  // [助詞, 「〜すると」の形]
  kmno4: ['を', '硫酸酸性の過マンガン酸カリウム水溶液で十分に酸化すると'],
  mild_oxidation: ['を', '硫酸酸性の二クロム酸カリウム水溶液で穏やかに酸化すると'],
  hydrolysis: ['を', '加水分解すると'],
  ozonolysis: ['を', 'オゾン分解すると'],
  kmno4_cleave: ['を', '硫酸酸性の過マンガン酸カリウム水溶液で酸化すると'],
  dehydration: ['に', '濃硫酸を加えて加熱し、分子内で脱水させると'],
  dehydration_ozonolysis: ['を', '脱水して得られるアルケンをオゾン分解すると'],
  hydrogenation: ['に', '白金触媒を用いて水素を十分に付加させると'],
  periodate: ['に', '過ヨウ素酸を十分に作用させると'],
  markovnikov: ['に', '酸触媒を用いて水を付加させると'],
  nitration: ['に', '濃硝酸と濃硫酸の混合物を作用させると'],
  bromination_fe: ['に', '鉄粉を触媒として臭素を作用させると'],
  chlorination: ['に', '鉄粉を触媒として塩素を作用させると'],
  sulfonation: ['に', '濃硫酸を加えて加熱すると'],
  bromine_water: ['に', '十分な量の臭素水を加えると'],
  nitro_reduction: ['を', 'スズと濃塩酸で還元し、水酸化ナトリウム水溶液で中和すると'],
  acetylation: ['に', '十分な量の無水酢酸を作用させると'],
  acetylation_primary: ['に', '同じ物質量の無水酢酸を作用させると'],
  deamination: ['を', '亜硝酸ナトリウムと塩酸で冷やしながらジアゾ化し、H₃PO₂ で還元すると'],
  diazo_hydrolysis: ['を', 'ジアゾ化し、その水溶液を温めると'],
  azo_coupling: ['を', 'ジアゾ化してナトリウムフェノキシド水溶液に加えると'],
  imide_hydrolysis: ['を', '穏やかな条件で加水分解すると'],
  ether_hydrogenolysis: ['を', '触媒の存在下で水素と反応させると'],
  ring_hydrogenolysis: ['を', '触媒を用いて水素と反応させると'],
  acetal_hydrolysis: ['に', '希酸を加えて加水分解すると'],
  acetal_etoh: ['を', '少量の硫酸を含む大過剰のエタノール中に十分な時間置くと'],
  acetal_meoh: ['を', '少量の硫酸を含む大過剰のメタノール中に十分な時間置くと'],
  acetonide: ['を', '酸触媒の存在下でアセトンと反応させると'],
  methylation_analysis: ['の', 'ヒドロキシ基をすべてメチル化してから加水分解すると'],
  bromine_addition: ['に', '臭素を付加させると'],
  br2_anti: ['に', '臭素を付加させると'],
  h2_syn: ['に', '白金触媒を用いて水素を付加させると'],
  nitric_oxidation: ['を', '硝酸で酸化すると'],
  sugar_degrade: ['から', '炭素を1つ減らす反応を行うと'],
};
// 「〜すると」→「〜しても」
function concessive(v) {
  return v.replace(/すると$/, 'しても').replace(/くと$/, 'いても').replace(/うと$/, 'っても').replace(/ると$/, 'ても');
}

function sentence(card, subject, result, cards) {
  if (S[card]) return S[card](subject, result);
  if (card === 'partial_hydrolysis') return `${subject} を途中まで加水分解すると、次の化合物が${result ? '得られた' : '得られなかった'}。`;
  if (card === 'acetylation_primary' && Array.isArray(result) && !result.length) return `${subject} は第一級アルコールの部分構造（–CH₂OH）をもたない。`;
  const [pa, verb] = V[card] || ['について', `「${cards[card].action}」の操作を行うと`];
  if (Array.isArray(result)) return result.length ? `${subject} ${pa}${verb}、次の化合物が得られた。` : `${subject} ${pa}${concessive(verb)}、変化はなかった。`;
  return `${subject} について「${cards[card].action}」を調べると、${result}。`;
}

// この結果からわかること（部分構造・位置関係）
const I = {
  silver_mirror: (r) => (r ? 'ホルミル基 –CHO をもつ（アルデヒド）。ギ酸エステル H–CO–O– も陽性' : '–CHO をもたない。C=O があるならケトン（鎖の途中の C=O）'),
  fehling: (r) => (r ? '–CHO をもつ（アルデヒド）' : '–CHO をもたない'),
  iodoform: (r) => (r ? 'CH₃–CO–C（メチルケトン）か CH₃–CH(OH)–C の部分構造をもつ。エタノールとアセトアルデヒドも陽性' : 'CH₃–CO– と CH₃–CH(OH)– のどちらももたない'),
  sodium: (r) => (r ? '–OH（アルコール・フェノール）か –COOH をもつ' : '–OH も –COOH ももたない（エーテル・エステル・ケトンなど）'),
  nahco3: (r) => (r ? '–COOH をもつ（炭酸より強い酸）。フェノールは溶けない' : '–COOH をもたない'),
  fecl3: (r) => (r ? 'ベンゼン環に直接 –OH がつく（フェノール類）' : 'フェノール性の –OH をもたない。ベンゼン環があっても OH は側鎖の炭素につくか、ない'),
  bromine: (r) => (r ? 'C=C か C≡C をもつ（ベンゼン環は臭素水と反応しない）' : 'C=C・C≡C をもたない。不飽和度は環・ベンゼン環・C=O による'),
  naoh: (r) => (r ? '酸性の基（–COOH かフェノール性 –OH）をもつ' : '酸性の基をもたない'),
  hcl: (r) => (r ? '塩基性の基（アミノ基）をもつ' : 'アミノ基をもたない（アミドの N は塩基性を示さない）'),
  chiral: (r) => (r ? `4つの異なる原子・原子団が結合した炭素が ${r} 個ある` : '4つの異なる原子・原子団が結合した炭素がない。対称な形か、枝の2本が同じ'),
  cis_trans: (r) => (r ? 'C=C の両端の炭素が、それぞれ異なる2つの基をもつ' : 'C=C がないか、C=C の一方の端に同じ基が2つつく（=CH₂、=C(CH₃)₂ など）'),
  anhydride: (r) => (r ? '2つの –COOH が近い（ベンゼン環のオルト位、C=C のシス側）' : '2つの –COOH が離れている（メタ・パラ位、トランス形）か、–COOH が1つ'),
  carbon_env: (r) => `分子を対称性で重ねると、炭素は ${r} 種類。等価な炭素の組を探すと骨格が絞れる（p-置換は対称性が高い）`,
  ring_cl: (r) => `ベンゼン環上の H の種類が ${r}。置換基の数と位置関係（オルト・メタ・パラ）が決まる`,
  cl_sub: (r) => `炭素に結合した H の種類が ${r}。鎖の枝分かれと対称性が決まる`,
  h2_uptake: (r) => `C=C が ${r} 個（C≡C なら1つで 2 mol）。ベンゼン環には付加しない`,
  dehydration_count: (r) => `OH の炭素の隣にある H をもつ炭素の数と、できる C=C のシス-トランスから、OH の位置がわかる（${r} 種類）`,
  stereo_count: (r) => `不斉炭素とシス-トランスの数から 2ⁿ 個が上限。メソ体があるとそれより少ない（${r} 種類）`,
  optically_active: (r) => (r ? '鏡像と重ならない（不斉炭素をもち、対称面がない）' : '鏡像と重なる（不斉炭素がないか、対称面をもつメソ体）'),
  kmno4: () => 'ベンゼン環に直結した炭素は、長さにかかわらず –COOH になる。生成物の COOH の数と位置 = 側鎖の数と位置。第一級アルコールはカルボン酸、第二級はケトンになる',
  mild_oxidation: () => '第一級アルコール → アルデヒド、第二級 → ケトン。酸化されなければ第三級アルコール',
  hydrolysis: () => 'エステル結合・アミド結合で切れた断片。H₂O の数 = 切れた結合の数',
  ozonolysis: () => 'C=C を切って両側を C=O にした断片。断片の C=O どうしをつなぎ直すと C=O の位置が C=C の位置',
  kmno4_cleave: () => 'C=C で切れ、H が2つの端は CO₂、H が1つの端は –COOH、H のない端はケトンになる',
  dehydration: () => 'OH の炭素と、その隣の H をもつ炭素の間に C=C ができる。生成物の C=C の位置から OH の位置がわかる',
  dehydration_ozonolysis: () => '脱水でできた C=C の位置で切れる。断片から OH の位置と骨格がわかる',
  hydrogenation: () => 'C=C が消えて炭素骨格（枝分かれ）が決まる。同じアルカンになる化合物は同じ骨格',
  periodate: () => '隣り合う OH（または C=O）の炭素の間の結合が切れる。生成物から OH の並びがわかる',
  markovnikov: () => 'OH は H の少ない側の炭素に、H は H の多い側の炭素につく',
  nitration: () => '–CH₃・–OH・–NH₂ などはオルト・パラ、–NO₂・–COOH・–CHO などはメタに入る。空いている位置の数も手がかり',
  bromination_fe: () => '配向性に従って1か所に入る。オルト・パラ配向かメタ配向かで置換基がわかる',
  chlorination: () => '配向性に従って1か所に入る',
  sulfonation: () => '配向性に従って1か所に入る',
  bromine_water: () => '–OH・–NH₂ のオルト位とパラ位の空いた場所に全部入る。入った Br の数 = 空いた場所の数',
  nitro_reduction: () => '–NO₂ が –NH₂ になる',
  acetylation: () => '–OH と –NH₂ が –OCOCH₃・–NHCOCH₃ になる。分子量が 42 ずつ増える',
  acetylation_primary: () => '第一級アルコール –CH₂OH だけがアセチル化される。どの OH が第一級かがわかる',
  deamination: () => 'ベンゼン環の –NH₂ が H に置き換わる',
  diazo_hydrolysis: () => 'ベンゼン環の –NH₂ が –OH に置き換わる',
  azo_coupling: () => 'フェノールの OH のパラ位（ふさがっていればオルト位）に –N=N– でつながる',
  imide_hydrolysis: () => 'C(=O)–N–C(=O) の C–N が1本だけ切れる。どちらが切れても生成物に含まれる',
  ether_hydrogenolysis: () => 'Ar–O–Ar′ の C–O が切れて、フェノールと芳香族炭化水素になる',
  ring_hydrogenolysis: () => '三員環・四員環の C–C が1本切れる。ひずみの残らない結合が切れる',
  acetal_hydrolysis: () => 'アセタールがカルボニル化合物とアルコールに戻る',
  acetal_etoh: () => '五員環・六員環をつくれるアルデヒドは環状アセタール、つくれないものはジエチルアセタールになる',
  acetal_meoh: () => '五員環・六員環をつくれるアルデヒドは環状アセタール、つくれないものはジメチルアセタールになる',
  acetonide: () => '近い2つの OH（1,2 位か 1,3 位）がアセトンと環をつくる。どの OH が近いかがわかる',
  methylation_analysis: () => 'メチル化されずに残った OH の位置 = グリコシド結合に使われていた位置',
  bromine_addition: () => 'C=C の両端に Br が1つずつつく',
  br2_anti: () => '2つの Br は C=C の平面の反対側から付加する。E と Z で生成物の立体が違う',
  h2_syn: () => '2つの H は平面の同じ側から付加する',
  nitric_oxidation: () => '両端が –COOH になる。対称面ができればメソ体（光学不活性）',
  sugar_degrade: () => 'C1 が外れ、C2 が –CHO になる。残りの不斉炭素の配置は変わらない',
  partial_hydrolysis: (r) => (r ? 'この断片がこのつながり方でエステル結合している' : 'このつながり方はしていない'),
};
function infer(card, result, cards) {
  if (I[card]) return I[card](result);
  return cards && cards[card] ? cards[card].action : '';
}

// 反応名（空欄補充の問い）
const NAMES = {
  silver_mirror: '銀鏡反応', fehling: 'フェーリング反応', iodoform: 'ヨードホルム反応', fecl3: '塩化鉄(III) による呈色反応',
  ninhydrin: 'ニンヒドリン反応', xanthoprotein: 'キサントプロテイン反応', biuret: 'ビウレット反応', ozonolysis: 'オゾン分解',
  nitration: 'ニトロ化', sulfonation: 'スルホン化', hydrolysis: '加水分解', acetylation: 'アセチル化', azo_coupling: 'ジアゾカップリング',
};

module.exports = { sentence, infer, NAMES };
