# Нестационарная ландауэровская эффективность — открытые материалы

*([English version](./README.md))*

Открытые материалы к статье **«The Non-Stationary Landauer Efficiency: Memory Growth, Informational Nostalgia, and the Thermodynamics of Learning Systems»** (русский — основной язык; включён английский перевод).

Это продолжение статьи *«Landauer Efficiency of Self-Modeling: An Operational Scale of Vitality»* ([andriishin/landauer-self-modeling-oa](https://github.com/andriishin/landauer-self-modeling-oa)), расширяющее ландауэровскую эффективность

```
η_L = I_pred / N_max,   N_max = E_actual / (k_B T ln 2)
```

со стационарного режима на **нестационарную** среду с растущей памятью. Система удерживает свою модель против среды, которая дрейфует быстрее, чем модель обновляется, поэтому растущая доля хранимых битов перестаёт что-либо предсказывать — *информационная ностальгия*. Доказаны два результата: **Лемма 1** (удержание бита против термодинамической эрозии требует строго положительной мощности с явной нижней границей) и **Лемма 2** (при медленном дрейфе Орнштейна–Уленбека и любом конечном темпе обновления доля ностальгии асимптотически остаётся выше явной положительной константы). Следствия: эффективность обнуляется выше критической ностальгии и существует оптимальный темп полного сброса памяти.

## Состав

```
paper/
  main.md / main.ru.md               Рукопись — английский (язык подачи) / русский (основной)
  supplementary.md / supplementary.ru.md   Дополнительные материалы — англ. / рус.
  main.tex / main.ru.tex             Исходник LaTeX (шаблон Springer Nature sn-jnl; генерируется из .md)
  supplementary.tex / supplementary.ru.tex
  main.pdf (48 с.) / main.ru.pdf (50 с.)            Скомпилированная рукопись
  supplementary.pdf (43 с.) / supplementary.ru.pdf (45 с.)
  refs.bib                           Библиография
  sn-jnl.cls, sn-*.bst               Шаблон Springer Nature + стили библиографии (third-party; см. LICENSE)
  figs/                              Рисунки (генерируются симуляциями ниже)
simulations/                         Воспроизводимые численные эксперименты
  markov_drift/                      PSP-суррогат (§ 6.1): три режима η_L(t)
  markov_drift_ou/                   OU-дрейф (§ 6.2): непрерывный дрейф логитов, режим коллапса
  markov_drift_ou_iinf/              OU + кумулятивная BNT-избыточность (логарифмический коэффициент joint-fit)
  markov_drift_ou_iinf_adiab/        Адиабатический скан: коэффициент BNT vs σ²/λ (Supplementary § S8.3)
  permuted_mnist/                    Контролируемый ML-камертон (§ 8.4): катастрофическое забывание
```

## Пересборка PDF

Готовые PDF включены для прямого чтения. Для пересборки из `.tex` нужен дистрибутив TeX с `pdflatex` и `bibtex` (TeX Live или MiKTeX). Файлы шаблона Springer Nature (`sn-jnl.cls`, `sn-*.bst`) вложены, ничего скачивать не требуется. Стиль ссылок — `sn-mathphys-ay` (автор–год). Стандартный цикл для любого документа, например `main`:

```bash
cd paper
pdflatex main.tex
bibtex   main
pdflatex main.tex
pdflatex main.tex
```

(замените `main` на `supplementary`, `main.ru` или `supplementary.ru` для остальных документов). Рисунки в `paper/figs/` нужны для компиляции и пересоздаются симуляциями.

## Воспроизведение симуляций

Каждый `simulations/<эксперимент>/` самодостаточен:

```bash
cd simulations/markov_drift_ou && pip install -r requirements.txt && python main.py
```

Инициализация зерном `SEED=20260524` с посерийными сдвигами; на той же версии NumPy каждый запуск воспроизводит рисунки и числа в `results_summary.{txt,json}`. В каждой папке есть `README.md` (английский) и `README.ru.md` (русский) с описанием модели, параметров и устанавливаемого результата.

## Лицензия

- **Текст, рисунки, данные** (рукопись, дополнительные материалы, библиография, рисунки, README симуляций, ожидаемые выводы): **Creative Commons Attribution 4.0 International (CC BY 4.0)** — см. [`LICENSE`](./LICENSE).
- **Исходный код** в `simulations/` (файлы `*.py`): **лицензия MIT** — см. [`simulations/LICENSE`](./simulations/LICENSE).

## Цитирование

Архивные материалы (этот репозиторий):

> Andriishin, A. *The Non-Stationary Landauer Efficiency: Memory Growth, Informational Nostalgia, and the Thermodynamics of Learning Systems.* Zenodo, 2026. https://doi.org/10.5281/zenodo.20653051

Депозит опубликован на Zenodo с указанным выше DOI. Машиночитаемое описание цитирования — в [`CITATION.cff`](./CITATION.cff). Журнальная версия записи будет добавлена отдельной ссылкой после принятия.
