# jouzoukettei — 有機構造決定パズル

京大化学の有機（大問3の構造決定）で使う「条件から候補を絞り込む推論」を反復訓練するブラウザゲーム。
仕様は Claude Docs の「有機構造決定パズル 仕様書」。

| モード | 内容 | 問題数 |
| --- | --- | --- |
| 絞り込み型 | 分子式と手がかりカードから、候補を消して1つに決める（C₅以下の脂肪族） | 20 |
| 大問型 | C₁₂〜C₂₀ のエステル・アミド・ペプチドを加水分解 → 切れた結合の数 → 断片ごとに絞り込み → 部分加水分解で組み立て | 8 |
| 数え上げ型 | 条件を満たす構造異性体と立体異性体の数を答え、全異性体で答え合わせ | 10 |
| 高分子型 | 平均分子量から重合度と結合の数を計算 → 加水分解で得た単量体を絞り込む（ナイロン66・ナイロン6・PET・ポリ乳酸・ポリ酢酸ビニル・アラミド） | 6 |

## 遊び方

`dist/index.html` をブラウザで開くだけ。サーバーもネット接続もいらない（フォントだけ Google Fonts から読む）。
`dist/` は GitHub Actions が自動で生成してコミットする。

1. 問題を選ぶ（難易度 ★1〜★3、または「おまかせ」）
2. 手がかりカードをタップして結果を見る。1枚ごとに −10 点
3. 矛盾する候補を自分で消す（左右スワイプか「消す」）。自動では消えない
4. 残り1つで「確定」。誤答なら矛盾したカード、根拠なしで消した候補は消すべきだったカードが出る
5. 間違えたカードの種類が「復習リスト」に溜まり、おまかせ出題で多めに出る

23時〜5時はプレイできない（仕様書の運用ルール）。

## 構成

| パス | 役割 |
| --- | --- |
| `problems/narrow.json` | 絞り込み型（仕様書の JSON 形式＋`level`） |
| `problems/big.json` | 大問型（`fragments` と `assemble`） |
| `problems/polymer.json` | 高分子型（`unit` は * で端を示した繰り返し単位、`n` は重合度、`bondsPerUnit` は縮合で外れる H₂O の数） |
| `problems/count.json` | 数え上げ型（`pool` はその分子式・分類の異性体すべて、`expect` は教科書の数） |
| `problems/names.json` | 結果画面に出す物質名 |
| `src/chem.js` | 手がかりカードの判定ロジック（SMARTS・酸化/加水分解/オゾン分解・不斉炭素・シス-トランス） |
| `scripts/test-rules.js` | 判定ロジックの単体テスト（教科書の典型例） |
| `scripts/validate.js` | 問題の自動検証 |
| `scripts/build.js` | 検証 → 構造式SVGと判定表を前計算 → `dist/index.html` を出力 |
| `src/index.html` | ゲーム画面のテンプレート |

化学処理（RDKit）はビルド時に Node で行い、ブラウザには前計算した構造式SVGと判定結果だけを埋め込む。
仕様書ではブラウザで RDKit.js を CDN から読む案だったが、wasm 読み込みが不要になり、スマホでも一瞬で開けてオフラインでも動くのでこちらにした。

## 問題を足す

```sh
npm install
npm test          # 判定ロジックのテスト + 全問題の検証
npm run build     # dist/index.html を作り直す
```

検証スクリプトが確認すること（1つでも落ちた問題は収録しない）:

1. 全候補の分子式が `formula` と一致する
2. 候補に重複がない（canonical SMILES で比較）
3. 全手がかりを適用すると、候補が `answer` の1つだけ残る
4. どの手がかりも、単独で候補を1つ以上消す
5. 手がかりの `result` が、正解の構造に実際にカードを当てた結果と一致する
6. カードが適用できない構造（例: C=C を持つ候補に KMnO₄）を含まない
7. 大問型: X の加水分解生成物と断片が一致し、原子の収支（X ＋ n H₂O）が合う。組み立ての候補はどれも同じ断片に分かれる
9. 高分子型: 繰り返し単位 ＋ `bondsPerUnit` H₂O と単量体の合計で原子の収支が合う
8. 数え上げ型: 母集団が重複なく所定の数そろっていて、答えの数が `expect`（教科書の数）と一致する

## 手がかりカード

| card | 画面の名前 | result |
| --- | --- | --- |
| `silver_mirror` | 銀鏡反応 | true / false |
| `fehling` | フェーリング液 | true / false |
| `iodoform` | ヨードホルム反応 | true / false |
| `sodium` | 金属Na | true / false |
| `nahco3` | NaHCO₃ | true / false |
| `fecl3` | FeCl₃ | true / false |
| `bromine` | 臭素水 | true / false |
| `kmno4` | KMnO₄ 酸化 | 生成物 SMILES の配列（`[]` = 酸化されない） |
| `hydrolysis` | 加水分解 | 生成物 SMILES の配列（`[]` = 加水分解されない） |
| `ozonolysis` | オゾン分解 | 生成物 SMILES の配列（`[]` = 反応しない） |
| `chiral` | 不斉炭素 | 個数 |
| `cis_trans` | シス-トランス異性 | true / false |
| `naoh` / `hcl` | NaOH 水溶液 / 希塩酸に溶けるか | true / false |
| `dehydration` | 分子内脱水の生成物 | SMILES の配列 |
| `dehydration_count` | 脱水で生じるアルケンの種類（シス・トランスを別に数える） | 数 |
| `dehydration_ozonolysis` | 脱水 → オゾン分解の生成物 | SMILES の配列 |
| `ring_cl` / `cl_sub` | ベンゼン環 / 炭素上の H を1つ Cl にした化合物の種類 | 数 |
| `anhydride` | 加熱で環状の酸無水物になるか | true / false |
| `partial_hydrolysis` | 部分加水分解でその化合物が得られるか | SMILES（`pep:Ala-Phe` 表記も可） |
| `ninhydrin` / `alpha_amino` / `xanthoprotein` / `sulfur` / `biuret` | アミノ酸・ペプチドの検出反応 | true / false |

ゲーム画面にはカードの結果しか出さない。「その結果から何が確定するか」は仕様書の手がかり表の右端の列に自分で書く。
