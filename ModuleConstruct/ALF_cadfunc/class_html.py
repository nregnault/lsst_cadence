#Def class for build html

class MakeHtml():

	def __init__(self, title, template={}):

		self.title = title

		if template == {}:
			template['init'] = '/home/alamure/cadence/software/lsst_cadence/ModuleConstruct/ALF_cadfunc/template_initial.html'
			template['head'] = '/home/alamure/cadence/software/lsst_cadence/ModuleConstruct/ALF_cadfunc/template_head.html'
			template['footer'] = '/home/alamure/cadence/software/lsst_cadence/ModuleConstruct/ALF_cadfunc/template_footer.html'
			template['ending'] = '/home/alamure/cadence/software/lsst_cadence/ModuleConstruct/ALF_cadfunc/template_ending.html'

		self.init = open(template['init'], 'r').read()
		self.head = open(template['head'], 'r').read()
		self.footer = open(template['footer'], 'r').read()
		self.ending = open(template['ending'], 'r').read()

		self.html = self.init + self.head

	def add_Hover(self, text, link, hover='\n'):

		hover += '<div class="w3-top"><div class="w3-bar w3-theme w3-top w3-left-align w3-large">'
		hover += '<a href="javascript:void(0)" class="w3-bar-item w3-button w3-right w3-hide-large w3-hover-white w3-large w3-theme-l1" onclick="w3_open()">'
		hover += '<i class="fa fa-bars"></i></a>'
		for text_i, link_i in zip(text, link):
			hover += '<a href="{}" class="w3-bar-item w3-button w3-hide-small w3-hover-white">{}</a>'.format(link_i, text_i)
		hover += '</div></div>\n'

		self.html += hover

	def add_Navigation(self, text, link, menu='Menu', nav='\n'):

		nav += '<nav class="w3-sidebar w3-bar-block w3-collapse w3-large w3-theme-l5 w3-animate-left" id="mySidebar">'
		nav += '<a href="javascript:void(0)" class="w3-right w3-xlarge w3-padding-large w3-hover-black w3-hide-large" onclick="w3_close()" title="Close Menu">'
		nav += '<i class="fa fa-remove"></i></a>'
		nav += '<h4 class="w3-bar-item"><b>{}</b></h4>\n'.format(menu)

		for text_i, link_i in zip(text, link):

			nav += '<a href="{}" class="w3-bar-item w3-button w3-hover-black">{}</a>\n'.format(link_i, text_i)

		nav += '<div style="cursor:pointer" title="close side menu" id="myOverlay" onclick="w3_close()" class="w3-overlay w3-hide-large"></div></nav>\n'

		self.html += nav

	def begin_body(self, body='\n'):

		body += '<div style="margin-left:250px" class="w3-main">'
		body += '<div class="w3-row w3-padding-64">'
		body += '<div class="w3-container">\n'

		self.html += body

	def ending_body(self, body='\n'):

		body += '</div></div></div>\n'

		self.html += body

	def add_label(self, label, h='h1'):

		self.html += '\n<{} style="text-align:center">{}</{}>\n'.format(h, label, h)

	def add_pict(self, link, pict='\n', width=1200):

		pict += '<div class="w3-center">'
		pict += '<img src="{}" style="float:center" width="{}px">'.format(link, str(width))
		pict += '</div>\n'

		self.html += pict

	def add_mp4(self, link, video='\n', width=80):

		video += '<div style="text-align:center">'
		video += '<video width="{}%" controls="1">'.format(str(width))
		video += '<source src="{}">'.format(link)
		video += '</video></div>'

		self.html += video

	def close(self):

		self.ending_body()
		self.html += self.footer
		self.html += self.ending

		f = open(self.title + '.html', 'w')
		f.write(self.html)
		f.close()
