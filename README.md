# Numerical Recipes applied to Thermodynamics site files

Jekyll and a mix of HMTL5 Boilerplate and my own stuff.

* prefix-free
* prism (for the line numbers to be drawn, the closing tags must be in their own line)
* MathJax
* Several fonts but optimized with `&text=NumericalRps` and `&text=THERMODYNAICS` to download just what I need. I write this here to remember to change it if I change the `text-transform` property.
* When copying code, change not only the &lt; but the ampersand too.
* submenus in the includes

**TODO:**
 - [ ] Check [http://couscous.io](http://couscous.io) out and see if it is better and if I can use it on GitHub.
 - [ ] Find a way to strip the `.html` from page URLs (see below).


### What I wanted for this site

Sometimes you want to sort stuff in a custom order. Jekyll sorts alphabetically, by date, etc. Also, it doesn't provide a way to have subcategories or even tags.
I guess Jekyll is just for very, very simple blogs.
Check this [SO link](http://stackoverflow.com/questions/27191110/frontmatter-automation-and-category-sorting-in-jekyll)


### Settings to achieve that

* *Folder structure*: Folders named after categories. Files inside them are pages under that category. Category folders under root. Put an `index.html` file under each folder to avoid 403 errors.

* *File naming*: `title.html`

* *Categories*: In the `_congif.yml` file, add:

```yml
categories:
  "breitwigner": "Breit-Wigner"
  "pvdiagrams":  "PV-Diagrams"
  "gibbs":       "Gibbs"
```

* *Listing*: Add this to the `nav` tag's list:

```liquid
{% for cat in site.categories %}
	<li><a href="{{ site.baseurl }}/{{ cat[0] }}/">{{ cat[1] }}</a></li>
{% endfor %}
```
* *Generated URLs*: `siteurl/category/title.html`. Looks like adding `permalink: /:categories/:title` to the `_congif.yml` file only works for posts but not for pages. I can not strip the `.html` from the page URL.

 
### Automated submenu generation

I also have some kind of "subpages" inside the main categories.

```html
{% assign sorted_pages = site.pages | sort:"date" %}
			<ul class="submenu">{% for cat in site.categories %}{% if page.submenu == cat[0] %}{% for post in sorted_pages %}{% if post.submenu == page.submenu %}
				<li{% if page.url == post.url %} class="active"{% endif %}><a href="{{ post.url }}">{{ post.stitle }}</a></li>{% endif %}{% endfor %}{% endif %}{% endfor %}
			</ul>
```

For this I need two things in the Front Matter:

* A `date` variable to sort them by date, otherwise Jekyll will sort them alphabetically.

* A short title variable (`stitle`) to place in the menu, since the actual title is too large to fit.

I sort the pages by date. Then I loop in the categories, if the category matches the directory of the current page, I enter a loop over all the pages and take only the ones that have the same directory or category.

In short:

```liquid
{% assign sorted_pages = site.pages | sort:"date" %}

{% for cat in site.categories %}
	{% if page.dir contains cat[0] %}
		{% for sorted in sorted_pages %}
			{% if sorted.dir == page.dir %}
				<li{% if page.url == sorted.url %} class="active"{% endif %}>
					<a href="{{ sorted.url }}">{{ sorted.stitle }}</a>
				</li>
			{% endif %}
		{% endfor %}
	{% endif %}
{% endfor %}
```

## License

MIT license.
